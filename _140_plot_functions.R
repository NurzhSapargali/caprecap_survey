# File: _140_plot_functions.R
#---------------------------#
# This script is contains functions to plot the results of capture-recapture simulations.
# The functions read in the simulation results, preprocess the data, aggregate
# it, and then generate plots to visualize the performance of different
# estimation methods under various conditions. The plots include metrics such as
# root mean squared error relative to the baseline, relative absolute bias,
# and log-transformed versions of these metrics

# Author: Nurzhan Sapargali
# Date created: 2024-02-26
#---------------------------#
library(tidyverse)
library(directlabels)

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
# - Gets the baseline values for "MPLE-NB"
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
  
  # Get the baseline values for "MPLE-NB"
  baseline <- agg[
    agg$type == "MPLE-NB",
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
  
  # Mark incomplete estimates
  agg$incomplete <- (agg$trials < TRIALS)
  incomp_by_pop <- agg[agg$incomplete,][c("type", "N")] %>%
    distinct(type, N)
  
  for (i in 1:nrow(incomp_by_pop)){
    estimator <- incomp_by_pop$type[i]
    pop <- incomp_by_pop$N[i]
    agg$type[agg$type == estimator & agg$N == pop] <- paste0(estimator, "*")
  }
  
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
  
  # Vector of types to exclude
  exclude_types <- c("MPLE-NB", "Morgan-Ridout", "Morgan-Ridout*")
  
  # Create the plot
  ggplot(
    # Get specified population and exclude "MPLE-NB" and "Morgan Ridout"
    agg[(agg$N == pop) & !(agg$type %in% exclude_types), ],
    mapping = aes(
      x = T,
      y = !!y_col,
      colour = type,
      group = type
    )
  ) +
    geom_line(alpha = 0.9) +
    geom_point(
      mapping = aes(colour = ifelse(incomplete, NA, type), shape = incomplete),
      size = 2.0
    ) +
    # Add a horizontal line at y = 0
    geom_hline(yintercept = 0, colour = "black", alpha = 0.35) +
    theme_minimal() +
    # Set the coordinate limits
    coord_cartesian(ylim = ylim, xlim = xlim) +
    # Set the decimal places for the y-axis
    scale_y_continuous(
      labels = scales::number_format(accuracy = 0.1),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    # Remove the colour guide
    guides(colour = "none", shape = "none") +
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


# 'comparison_plots' function for creating line plots of the aggregated data
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

comparison_plots <- function(
  agg,
  pop,
  ylim_bias = c(NA, NA),
  ylim_rmse = c(NA, NA),
  xlim = c(NA, 12.2),
  dense = NA
) {
  # Create the plot title based on the population and graph density
  title <- paste0("N = ", as.character(pop))
  
  # Create the RMSE plot
  p1 <- plot_lines(agg, pop, "lrel_rmse", ylim_rmse, "LogRelRMSE", title, xlim)
  # Create the bias plot
  p2 <- plot_lines(agg, pop, "lrel_bias", ylim_bias, "LogRelAbsBias", "", xlim)
  
  # Return the plots as a list
  list(p1, p2)
}


# 'estimates_boxplot' function for creating a box plot of the estimates
#
# This function creates a box plot of the estimates 'N_hat' for a specific
# population 'pop'. The true value of the population is specified by 'true_val'.
# The y-axis is labeled with 'ylab', and the variable to be plotted on the
# y-axis is specified by 'to_plot', with a default value of "N_hat".
#
# Args:
#   pr_results: A data frame containing the results to be plotted.
#   pop: The population to be plotted.
#   true_val: The true value of the population.
#   ylab: The label for the y-axis.
#   to_plot: The variable to be plotted on the y-axis. Default is "N_hat".
# Returns:
#   A ggplot object representing the box plot.
estimates_boxplot <- function(
  pr_results,
  pop,
  ylab,
  to_plot = "N_hat"
) {
  
  dummy_df <- pr_results[1:2,]
  
  pr_results %>%
    filter(N == pop) %>%
    ggplot(mapping = aes(x = T, y = !!sym(to_plot), fill = type)) +
    geom_boxplot() +
    stat_summary(
      fun = mean,
      color = "gold",
      position = position_dodge(0.75),
      geom = "point",
      shape = 17,
      size = 2.75,
      show.legend = FALSE,
    ) +
    geom_point(
      data = dummy_df,  # dummy data
      aes(shape = "Mean"),
      color = "gold",
      size = 2.75,
      alpha = 0.0
    ) +
    scale_shape_manual(
      values = c("Mean" = 17)
    ) +
    geom_hline(aes(yintercept = log(N), linetype = "True log N")) +
    scale_linetype_manual(
      name = "",
      values = c("True log N" = "dashed"),
      labels = c("True log N")
    ) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    theme_minimal() +
    theme_pubr() +
    guides(
      fill = "legend",
      linetype = "legend",
      shape = guide_legend(
        override.aes = list(alpha = 1, color = "gold", size = 2.75)
      )
    ) +
    theme(legend.title = element_blank()) +
    xlab("T") +
    ylab(ylab)
}