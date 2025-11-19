# File: _140_plot_functions.R
#---------------------------#
# This script contains functions to plot the results of simulations.
# The functions read in the simulation results, preprocess the data, aggregate
# it, and then generate plots to visualize the performance of different
# estimation methods under various conditions. The plots include metrics such as
# root mean squared error relative to the baseline, relative absolute bias,
# and log-transformed versions of these metrics

# Author: Nurzhan Sapargali
# Date created: 2024-02-26
#---------------------------#
library(dplyr)
library(ggplot2)
library(ggpubr)
library(directlabels)

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
  "Turing Geometric",
  "MPLE-G"
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
  "TG",
  "MPLE-G"
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
#   new_names: A vector of new names to replace the old names in the 'type'
# column.
#   old_names: A vector of old names to be replaced in the 'type'
#   minT: The minimum value of 'T' to include in the results. Rows with 'T' less
# than 'minT' are filtered out.
#
# Returns:
#   The preprocessed 'results' data frame.

preprocess <- function(results, new_names, old_names, minT = NA) {
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
  name_map <- setNames(new_names, old_names)

  # Use the named vector to replace the old names with the new names in the
  # 'type' column
  names_to_replace <- results[results$type %in% old_names, c("type")]
  results[results$type %in% old_names, c("type")] <- name_map[names_to_replace]

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
  agg <- pr_results |>
    group_by(type, T, N) |>
    summarise(
      mse = mean(sq_dev),
      hat_mean = mean(N_hat),
      trials = length(N_hat),
      hat_var = var(N_hat),
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
  ) |>
    mutate(
      lrel_bias = log(rel_bias),
      lrel_rmse = log(rel_rmse),
      lrel_var = log(rel_var)
    )

  # Remove the baseline columns
  agg <- select(agg, -ends_with("_baseline"))

  agg$incomplete <- FALSE
  for (t in levels(agg$T)) {
    for (n in unique(agg$N)) {
      max_trials <- max(
        agg$trials[agg$T == t & agg$N == n]
      )
      agg$incomplete[agg$T == t & agg$N == n & agg$trials < max_trials] <- TRUE
    }
  }
  # Mark incomplete estimates
  incomp_by_pop <- agg[agg$incomplete,][c("type", "N")] |>
    distinct(type, N)

  for (i in seq_len(nrow(incomp_by_pop))) {
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
    xlab("Sample draws") +
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
# by 'xlim', with a default range of NA to 15.0.
#
# Args:
#   agg: A data frame containing the aggregated data to be plotted.
#   pop: The population to be plotted.
#   ylim_bias: The limits for the y-axis of the bias plot.
# Default is c(-3.0, 3.0).
#   ylim_rmse: The limits for the y-axis of the RMSE plot.
# Default is c(-3.0, 3.0).
#   xlim: The limits for the x-axis. Default is c(NA, 15.0).
#
# Returns:
#   A list containing two ggplot objects representing the line plots.

comparison_plots <- function(
  agg,
  pop,
  ylim_bias = c(NA, NA),
  ylim_rmse = c(NA, NA),
  xlim = c(NA, 12.2)
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
# The plot includes a horizontal line indicating the true log N value.
# The mean of the estimates is highlighted with a gold triangle and a legend is
# included to indicate this by creating a dummy data point with alpha = 0.
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

  pr_results |>
    filter(N == pop) |>
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
    xlab("Sample draws") +
    ylab(ylab)
}


# 'rel_bias_plots' function for creating relative bias plots

# This function creates relative bias plots from the provided bias data
# and saves the plot to the specified output filename.
#
# Args:
#   bias_data: A list of data frames containing bias data.
#   output_filename: The filename where the plot will be saved.
# Returns:
#   A ggplot object representing the relative bias plot.

rel_bias_plots <- function(bias_data, output_filename) {
  bias_data <- do.call(rbind, bias_data) |>
    arrange(desc(type))

  bias_data$rel_bias <- (bias_data$hat_mean - bias_data$N) / bias_data$N
  bias_data$N <- as.factor(bias_data$N)

  ggplot(
    data = bias_data,
    mapping = aes(
      x = T,
      y = rel_bias,
      color = N,
      linetype = type,
      group = interaction(N, type)
    )
  ) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(
      x = "Sample draws",
      y = "Relative Bias",
      linetype = NULL
    ) +
    facet_wrap(
      vars(alpha),
      labeller = labeller(alpha = function(x) paste("alpha =", x))
    ) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_linetype_manual(
      values = c("MPLE-NB" = "solid", "MPLE-G" = "dashed")
    ) +
    theme(legend.position = "top", legend.title = element_blank()) +
    theme_pubr()
}

# 'plot_everything' function for plotting all results
# This function generates and saves various plots based on the simulation
# results for different alpha values and population sizes.
# Args:
#   alphas: A vector of alpha values to process.
#   pops: A vector of population sizes to process.
#   new_names: A vector of new names for the estimation methods.
#   old_names: A vector of old names for the estimation methods.
#   estimates_folder: The folder path where the estimates CSV files are stored.
#   figures_folder: The folder path where the generated figures will be saved.
#   minT: The minimum value of 'T' to include in the results. Default is NA.
#   filename_suffix: A suffix to append to the filenames of the saved plots.
#   box_facet_size: A vector specifying the size of the boxplot facet plot
# (height, width) in mm
#   comp_facet_size: A vector specifying the size of the RMSE comparison facet
# plot (height, width) in mm
#   rel_bias_facet_size: A vector specifying the size of the relative bias facet
# plot (height, width) in mm
#
# Returns:
#   None. The function saves the generated plots to the specified folder.

plot_everything <- function(
  alphas,
  pops,
  new_names,
  old_names,
  estimates_folder,
  figures_folder,
  minT = NA,
  filename_suffix = "",
  box_facet_size = c(297, 210),
  comp_facet_size = c(297, 210),
  rel_bias_facet_size = c(99, 210),
  ylim_bias = c(NA, NA),
  ylim_rmse = c(NA, NA),
  xlim_comp = c(NA, 12.2)
) {
  bias_data <- list()

  # RMSE and absolute bias and boxplots of estimates for each alpha and N
  for (het in alphas){
    het_str <- format(het, nsmall = 1)
    prep <- read.csv(
      paste0(estimates_folder, "estimates_", het_str, filename_suffix, ".csv")
    ) |>
      preprocess(new_names, old_names, minT = minT)

    agg <- aggregate_data(prep)
    bias_data[[het_str]] <- agg[
      agg$type %in% c("MPLE-NB", "MPLE-G"), c("N", "T", "hat_mean", "type")
    ]
    bias_data[[het_str]]$alpha <- het

    comp_ps <- c()
    box_ps <- list()
    prep$logN_hat <- log(prep$N_hat)

    for (i in seq_along(pops)) {
      N <- pops[i]
      N_str <- as.character(N)
      comp_ps <- c(
        comp_ps,
        comparison_plots(
          agg, N, ylim_bias = ylim_bias, ylim_rmse = ylim_rmse, xlim = xlim_comp
        )
      )

      box_plot <- estimates_boxplot(
        prep[prep$type %in% c("MPLE-G", "MPLE-NB"),],
        N,
        "Log population estimate",
        "logN_hat"
      ) + ggtitle(paste0("N = ", N))
      box_ps[[i]] <- box_plot
    }

    # Collate and save comparison plots
    figure_rows <- length(pops)
    ggarrange(plotlist = comp_ps, ncol = 2, nrow = figure_rows)
    ggsave(
      paste0(figures_folder, "estimates_", het_str, filename_suffix, ".pdf"),
      width = comp_facet_size[2],
      height = comp_facet_size[1],
      units = "mm"
    )

    # Collate and save box plots
    ggarrange(
      plotlist = box_ps,
      ncol = 1,
      nrow = figure_rows,
      common.legend = TRUE,
      legend = "top"
    )
    ggsave(
      paste0(
        figures_folder, "estimates_box_", het_str, filename_suffix, ".pdf"
      ),
      width = box_facet_size[2],
      height = box_facet_size[1],
      units = "mm"
    )
  }

  # Plot relative bias
  bias_data <- do.call(rbind, bias_data) |>
    arrange(desc(type))
  bias_data$rel_bias <- (bias_data$hat_mean - bias_data$N) / bias_data$N
  bias_data$N <- as.factor(bias_data$N)

  ggplot(
    data = bias_data,
    mapping = aes(
      x = T,
      y = rel_bias,
      color = N,
      linetype = type,
      group = interaction(N, type)
    )
  ) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    labs(
      x = "Sample draws",
      y = "Relative Bias",
      linetype = NULL
    ) +
    facet_wrap(
      vars(alpha),
      labeller = labeller(alpha = function(x) paste("alpha =", x))
    ) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_linetype_manual(
      values = c("MPLE-NB" = "solid", "MPLE-G" = "dashed")
    ) +
    theme(legend.position = "top", legend.title = element_blank()) +
    theme_pubr()

  ggsave(
    paste0(figures_folder, "relative_bias", filename_suffix, ".pdf"),
    width = rel_bias_facet_size[2],
    height = rel_bias_facet_size[1],
    units = "mm"
  )
}