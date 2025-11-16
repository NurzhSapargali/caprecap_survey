source("./_140_plot_functions.R")

library(ggpubr)

RESULTS_FOLDER <- "./_900_output/data/appendix/low_sample/" # Folder containing the simulation results
OUTPUT_FOLDER <- "./_900_output/figures/appendix/" # Folder to save the figures
TRIALS <- 1000 # Number of trials in the simulation
ALPHAS <- c(0.5, 2.0) # Different values of alpha in the simulation
POP_SIZES <- c(600, 1000, 5000) # Different population sizes in the simulation
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


# Plot the results for unequal probability capture-recapture
bias_data <- list()
for (het in ALPHAS){
  het_str <- format(het, nsmall = 1)
  prep <- read.csv(
    paste0(RESULTS_FOLDER, "estimates_", het_str, "_low_sample.csv")
  ) %>%
    preprocess()
  agg <- aggregate_data(prep)
  bias_data[[het_str]] <- agg[agg$type %in% c("MPLE-NB", "MPLE-G"), c("N", "T", "hat_mean", "type")]
  bias_data[[het_str]]$alpha <- het
  
  comp_ps <- c()
  box_ps <- list()
  prep$logN_hat <- log(prep$N_hat)
  
  for (i in seq_along(POP_SIZES)) {
    N <- POP_SIZES[i]
    N_str <- as.character(N)
    
    if ((N == 5000)&(het == 0.5)) {
      p <- comparison_plots(agg, N, ylim_bias = c(NA, 9.0))
    } else if ((N == 5000)&(het == 2.0)) {
      p <- comparison_plots(agg, N, ylim_bias = c(NA, 10.0))
    } else {
      p <- comparison_plots(agg, N)
    }
    comp_ps <- c(comp_ps, p)
    
    box_plot <- estimates_boxplot(
      prep[prep$type %in% c("MPLE-G", "MPLE-NB"),],
      N,
      "Log population estimate",
      "logN_hat"
    ) + ggtitle(paste0("N = ", N))
    box_ps[[i]] <- box_plot
  }
  
  ggarrange(plotlist = comp_ps, ncol = 2, nrow = 3)
  ggsave(
    paste0(OUTPUT_FOLDER, "estimates_", het_str, "_low_sample.pdf"),
    width = 210,
    height = 297,
    units = "mm"
  )
  
  ggarrange(plotlist = box_ps, ncol = 1, nrow = 3, common.legend = TRUE, legend = "top")
  ggsave(
    paste0(OUTPUT_FOLDER, "estimates_box_", het_str, "_low_sample.pdf"),
    width = 210,
    height = 297,
    units = "mm"
  )
}

bias_data <- do.call(rbind, bias_data) |>
  arrange(desc(type))
bias_data$rel_bias <- (bias_data$hat_mean - bias_data$N) / bias_data$N
bias_data$N <- as.factor(bias_data$N)


ggplot(
  data = bias_data,
  mapping = aes(x = T, y = rel_bias, color = N, linetype = type, group = interaction(N, type))
) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  labs(
    x = "Sample draws",
    y = "Relative Bias",
    linetype = NULL
  ) +
  facet_wrap(vars(alpha), labeller = labeller(alpha = function(x) paste("alpha =", x))) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_linetype_manual(values = c("MPLE-NB" = "solid", "MPLE-G" = "dashed")) +
  theme(legend.position = "top", legend.title = element_blank()) +
  theme_pubr()

ggsave(
  paste0(OUTPUT_FOLDER, "relative_bias_low_sample.pdf"),
  width = 210,
  height = 99,
  units = "mm"
)
