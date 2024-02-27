# File: 400_plot_simulation_results.R
#---------------------------#
# This script is designed to plot the results of capture-recapture simulations.
# The script reads in the simulation results, preprocesses the data, aggregates
# it, and then generates plots to visualize the performance of different
# estimation methods under various conditions. The plots include metrics such as
# root mean squared error relative to the baseline, relative absolute bias,
# and log-transformed versions of these metrics

# Author: Nurzhan Sapargali
# Date created: 2024-02-26
#---------------------------#
library(tidyverse)
library(directlabels)
library(ggpubr)

OUTPUT_FOLDER <- "./_900_output/figures/" # Folder to save the figures
OLD_NAMES <- c(
  "Chao Lee Jeng 0",
  "Chao Lee Jeng 1",
  "Chao Lee Jeng 2",
  "Conway-Maxwell-Poisson",
  "Jackknife k = 1",
  "Jackknife k = 2",
  "Jackknife k = 3",
  "Jackknife k = 4",
  "Jackknife k = 5",
  "Turing",
  "Turing Geometric"
) # Names of the methods in the simulation results
NEW_NAMES <- c(
  "SC,0",
  "SC,1",
  "SC,2",
  "LCMP",
  "JK,1",
  "JK,2",
  "JK,3",
  "JK,4",
  "JK,5",
  "TB",
  "TG"
) # New names for the methods in plots


# 'preprocess' function for preprocessing the results data frame
#
# This function performs several preprocessing steps on the 'results' data
# frame:
# - Filters out rows where 'T' is less than 'minT'
# - Filters out rows where 'N_hat' is NA, NaN, Inf, or less than 0
# - Adds a new column 'sq_dev' which is the square of the difference between
# 'N_hat' and 'N'.
# - Converts 'T' to a factor
# - Replaces the old names in the 'type' column with new names
#
# Args:
#   results: A data frame containing the results to be preprocessed.
#
#   minT: The minimum value of 'T' to include in the results. Rows with 'T' less
# than 'minT' are filtered out.
#
# Returns:
#   The preprocessed 'results' data frame.

preprocess <- function(results, minT = NA) {
  # If 'minT' is not NA, filter out rows where 'T' is less than 'minT'
  if (!is.na(minT)){
    results <- results[results$T >= minT, ]
  }

  # Filter out rows where 'N_hat' is NA, NaN, Inf, or less than 0
  results <- results[(!is.na(results$N_hat)) & (!is.nan(results$N_hat)), ]
  results <- results[(results$N_hat != Inf) & (results$N_hat >= 0), ]

  # Add a new column 'sq_dev' which is the square of the difference between
  #'N_hat' and 'N'
  results$sq_dev <- (results$N_hat - results$N)^2

  # Convert 'T' to a factor
  results$T <- as.factor(results$T)

  # Create a named vector to map the old names to the new names
  name_map <- setNames(NEW_NAMES, OLD_NAMES)

  # Use the named vector to replace the old names with the new names in the
  # 'type' column
  names_to_replace <- results[results$type %in% OLD_NAMES, c("type")]
  results[results$type %in% OLD_NAMES, c("type")] <- name_map[names_to_replace]

  # Return the preprocessed 'results' data frame
  results
}


# 'aggregate_data' function for aggregating the simulation results
#
# This function performs several steps to aggregate the 'pr_results' data frame:
# - Groups by 'type', 'T', and 'N' and calculates several summary statistics
# - Calculates the root mean square error (RMSE) and bias
# - Gets the baseline values for "Pseudolikelihood"
# - Joins the baseline values with the other methods
# - Calculates the relative metrics
# - Checks if all trials are complete
# - Removes the baseline columns
#
# Args:
#   pr_results: A data frame containing the simulation results to be aggregated.
#
# Returns:
#   The aggregated 'pr_results' data frame.

aggregate_data <- function(pr_results) {
  # Group by 'type', 'T', and 'N' and calculate several summary statistics
  agg <- pr_results %>%
    group_by(type, T, N) %>%
    summarise(
      mse = mean(sq_dev),
      hat_mean = mean(N_hat),
      trials = length(N_hat),
      hat_var = var(N_hat)
    )

  # Calculate the root mean square error (RMSE) and bias
  agg$rmse <- sqrt(agg$mse)
  agg$bias <- abs(agg$hat_mean - agg$N)

  # Get the baseline values for "Pseudolikelihood"
  baseline <- agg[
    agg$type == "Pseudolikelihood",
    c("T", "N", "rmse", "bias", "hat_var")
  ]

  # Join the baseline values with the other methods
  agg <- left_join(agg, baseline, by = c("T", "N"), suffix = c("", "_baseline"))

  # Calculate the relative metrics
  agg <- mutate(
    agg,
    rel_bias = bias / bias_baseline,
    rel_var = hat_var / hat_var_baseline,
    rel_rmse = rmse / rmse_baseline,
  ) %>%
    mutate(
      lrel_bias = log(rel_bias),
      lrel_rmse = log(rel_rmse),
      lrel_var = log(rel_var)
    )

  # Remove the baseline columns
  agg <- select(agg, -ends_with("_baseline"))

  # Return the aggregated 'pr_results' data frame
  agg
}


# 'plot_lines' function for creating a line plot of the aggregated data
#
# This function creates a line plot of the aggregated data 'agg' for a specific
# population 'pop'. The y-axis variable is specified by 'y', and the y-axis
# limits are specified by 'ylim'. The x-axis limits are specified by 'xlim',
# with a default range of NA to 15.0. The plot is labeled with 'ylab' on the
# y-axis and 'title' as the title.
#
# Args:
#   agg: A data frame containing the aggregated data to be plotted.
#   pop: The population to be plotted.
#   y: The variable to be plotted on the y-axis.
#   ylim: The limits for the y-axis.
#   ylab: The label for the y-axis.
#   title: The title for the plot.
#   xlim: The limits for the x-axis. Default is c(NA, 15.0).
#
# Returns:
#   A ggplot object representing the line plot.

plot_lines <- function(agg, pop, y, ylim, ylab, title, xlim) {
  # Convert the y variable to a symbol
  y_col <- sym(y)

  # Create the plot
  ggplot(
    # Get specified population and exclude "Pseudolikelihood"and "Morgan Ridout"
    agg[(agg$N == pop) & (agg$type != "Pseudolikelihood"), ],
    mapping = aes(
      x = T,
      y = !!y_col,
      colour = type,
      group = type
    )
  ) +
    geom_line() +
    geom_point() +
    # Add a horizontal line at y = 0
    geom_hline(yintercept = 0, colour = "black", alpha = 0.35) +
    theme_minimal() +
    # Set the coordinate limits
    coord_cartesian(ylim = ylim, xlim = xlim) +
    # Remove the colour guide
    guides(colour = "none") +
    xlab("T") +
    theme_pubr() +
    ylab(ylab) +
    # Set the plot title
    ggtitle(title) +
    # Add direct labels to the lines
    geom_dl(
      aes(label = type),
      method = list(
        dl.trans(x = x * 1.025, y = y * 1.0), 
        "last.points",
        "bumpup",
        cex = 0.8
      )
    )
}


# 'plot_results' function for creating line plots of the aggregated data
#
# This function creates two line plots of the aggregated data 'agg' for specific
# population 'pop'. The y-axis limits for the bias and RMSE plots are specified
# by 'ylim_bias' and 'ylim_rmse', respectively. The x-axis limits are specified
# by 'xlim', with a default range of NA to 15.0. The graph density is specified
# by 'dense', with a default value of NA.
#
# Args:
#   agg: A data frame containing the aggregated data to be plotted.
#   pop: The population to be plotted.
#   ylim_bias: The limits for the y-axis of the bias plot.
# Default is c(-3.0, 3.0).
#   ylim_rmse: The limits for the y-axis of the RMSE plot.
# Default is c(-3.0, 3.0).
#   xlim: The limits for the x-axis. Default is c(NA, 15.0).
#   dense: The graph density. Default is NA.
#
# Returns:
#   A list containing two ggplot objects representing the line plots.

plot_results <- function(
  agg,
  pop,
  ylim_bias = c(-3.0, 3.0),
  ylim_rmse = c(-3.0, 3.0),
  xlim = c(NA, 12.0),
  dense = NA
) {
  # Create the plot title based on the population and graph density
  if (is.na(dense)){
    title <- paste0("N = ", as.character(pop))
  } else {
    title <- paste0(
      "N = ", as.character(pop),
      ", graph density = ",
      as.character(dense)
    )
  }

  # Create the RMSE plot
  p1 <- plot_lines(agg, pop, "lrel_rmse", ylim_rmse, "LogRelRMSE", title, xlim)
  # Create the bias plot
  p2 <- plot_lines(agg, pop, "lrel_bias", ylim_bias, "LogRelAbsBias", "", xlim)

  # Return the plots as a list
  list(p1, p2)
}

# Plot the results for unequal probability capture-recapture
# Define the lookup table for y-axis limits
ylim_lookup <- list(
  "0.5" = list(
    "1000" = list("rmse" = c(NA, 1.3), "bias" = c(NA, 5.8)),
    "5000" = list("rmse" = c(NA, 2.0), "bias" = c(NA, 5.5)),
    "10000" = list("rmse" = c(NA, 2.0), "bias" = c(NA, 5.0))
  ),
  "1.0" = list(
    "1000" = list("rmse" = c(NA, 1.4), "bias" = c(NA, 3.0)),
    "5000" = list("rmse" = c(NA, 1.8), "bias" = c(NA, 5.8)),
    "10000" = list("rmse" = c(NA, 2.75), "bias" = c(NA, 6.5))
  ),
  "5.0" = list(
    "1000" = list("rmse" = c(NA, 3.0), "bias" = c(NA, 6.5)),
    "5000" = list("rmse" = c(NA, 1.25), "bias" = c(NA, 6.75)),
    "10000" = list("rmse" = c(NA, 1.25), "bias" = c(NA, 5.0))
  ),
  "10.0" = list(
    "1000" = list("rmse" = c(NA, 3.0), "bias" = c(NA, 10.0)),
    "5000" = list("rmse" = c(NA, 1.5), "bias" = c(NA, 7.5)),
    "10000" = list("rmse" = c(NA, 1.25), "bias" = c(NA, 7.5))
  ),
  "Inf" = list(
    "1000" = list("rmse" = c(NA, 2.75), "bias" = c(NA, 3.0)),
    "5000" = list("rmse" = c(NA, 1.5), "bias" = c(NA, 2.5)),
    "10000" = list("rmse" = c(NA, 1.25), "bias" = c(NA, 2.75))
  )
)

for (het in c("0.5", "1.0", "5.0", "10.0", "Inf")){
  agg <- read.csv(
    paste0("./_900_output/data/simulated/estimates_", het, ".csv")
  ) %>%
    preprocess() %>%
    aggregate_data()

  ps <- c()

  for (N in c(1000, 5000, 10000)) {
    ylim_bias <- ylim_lookup[[het]][[as.character(N)]][["bias"]]
    ylim_rmse <- ylim_lookup[[het]][[as.character(N)]][["rmse"]]

    ps <- c(
      ps,
      plot_results(agg, N, ylim_bias = ylim_bias, ylim_rmse = ylim_rmse)
    )
  }

  ggarrange(plotlist = ps, ncol = 2, nrow = 3)
  ggsave(
    paste0(OUTPUT_FOLDER, "simulated/estimates_", het, ".pdf"),
    width = 210,
    height = 297,
    units = "mm"
  )

}

# Plot the results for Barabasi-Albert graphs
# Define the lookup table for y-axis limits
ylim_lookup <- list(
  "0.025" = list("rmse" = c(NA, 2.0), "bias" = c(NA, 4.5)),
  "0.05" = list("rmse" = c(NA, 2.5), "bias" = c(NA, 4.0)),
  "0.1" = list("rmse" = c(NA, 1.75), "bias" = c(NA, 6.0))
)

ps <- c()
N <- 10000
for (dense in c(0.025, 0.05, 0.1)) {
  agg <- read.csv(
    paste0("./_900_output/data/graphs/estimates_ba_", dense, ".csv")
  ) %>%
    preprocess() %>%
    aggregate_data()

  ylim_bias <- ylim_lookup[[as.character(dense)]][["bias"]]
  ylim_rmse <- ylim_lookup[[as.character(dense)]][["rmse"]]

  ps <- c(
    ps,
    plot_results(
      agg,
      N,
      dense = dense,
      ylim_bias = ylim_bias,
      ylim_rmse = ylim_rmse
    )
  )
}

ggarrange(plotlist = ps, ncol = 2, nrow = 3)
ggsave(
  paste0(OUTPUT_FOLDER, "graphs/estimates_ba.pdf"),
  width = 210,
  height = 297,
  units = "mm"
)
