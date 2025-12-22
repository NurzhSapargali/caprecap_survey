library(tidyverse)
library(ggpubr)

SIM_FOLDER <- "./_900_output/data/simulated/"
APP_FOLDER <- "./_900_output/data/appendix/beta_bin/"
OUTPUT_FOLDER <- "./_900_output/figures/appendix/"
POP_SIZES <- c(1000, 5000, 10000)
ALPHAS <- c(0.5, 2.0)
INTERMEDIATE <- FALSE # Whether to consider file with intermediate results or final results

if (!dir.exists(OUTPUT_FOLDER)) {
  ok <- dir.create(OUTPUT_FOLDER, recursive = TRUE)
  if (!ok) stop("Failed to create directory: ", dir_path)
}

aggregate_data <- function(results) {
  # Group by 'T', and 'N' and calculate several summary statistics
  agg <- results %>%
    group_by(T, N) %>%
    summarise(
      hat_mean = mean(N_hat),
      hat_median = median(N_hat),
      hat_std = var(N_hat)^0.5
    )

  # Return the aggregated 'pr_results' data frame
  agg
}


for (het in ALPHAS){
  het_str <- format(het, nsmall = 1)
  
  filename <- paste0(APP_FOLDER, "estimates_", het_str, "_betabin.csv")
  if (INTERMEDIATE){
    filename <- paste0(APP_FOLDER, "estimates_", het_str, "_betabin_intermediate.csv")
  }
  
  app <- read.csv(filename)
  app$T <- as.factor(app$T)

  filename <- paste0(SIM_FOLDER, "estimates_", het_str, ".csv")
  if (INTERMEDIATE){
    filename <- paste0(SIM_FOLDER, "estimates_", het_str, "_intermediate.csv")
  }
  
  sim <- read.csv(filename)
  sim$T <- as.factor(sim$T)
  sim <- sim[sim$type == "MPLE-NB",]
  
  # If considering full results, match the simulation trial ids in both datasets
  if (!INTERMEDIATE){
    sim <- sim[sim$trial %in% app$trial,]
  }

  sim <- aggregate_data(sim)
  app <- aggregate_data(app)

  sim$type <- "MPLE-NB"
  app$type <- "Beta-binomial"

  prep <- rbind(sim, app)
  prep$N <- as.factor(prep$N)

  line_ps <- list()

  for (N in POP_SIZES) {
    N_str <- as.character(N)

    line_plot <- ggplot(
      prep[prep$N == N,],
      aes(x = T, y = hat_mean, group = type, color = type)
    ) +
      geom_line() +
      geom_point() +
      theme_minimal() +
      theme_pubr() +
      xlab("T") +
      ylab("Mean estimate") +
      geom_hline(yintercept = N, linetype = "dashed") +
      guides(color = "none") +
      ggtitle(paste0("N = ", N))

    if (N == 1000) {
      line_plot <- line_plot +
        guides(color = "legend") +
        theme(legend.title = element_blank(), legend.position = c(0.9, 0.9))
    }

    line_ps[[N_str]] <- line_plot
  }
  ggarrange(plotlist = line_ps, ncol = 1, nrow = 3)
  # Figures 2 and 3 in the supplementary material
  if (INTERMEDIATE){
    ggsave(
      paste0(OUTPUT_FOLDER, "instability_", het_str, "_intermediate.pdf"),
      width = 210,
      height = 297,
      units = "mm"
    )
  } else{
    ggsave(
      paste0(OUTPUT_FOLDER, "instability_", het_str, ".pdf"),
      width = 210,
      height = 297,
      units = "mm"
    )
  }
}
