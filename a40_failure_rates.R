library(tidyverse)
library(ggpubr)

RESULTS_FOLDER <- "./_900_output/data/simulated/" # Folder containing the simulation results
OUTPUT_FOLDER <- "./_900_output/figures/appendix/" # Folder to save the figures
TRIALS <- 1000 # Number of trials in the simulation
ALPHAS <- c(0.5, 2.0) # Different values of alpha in the simulation
POP_SIZES <- c(1000, 5000, 10000) # Different population sizes in the simulation
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

for (het in ALPHAS){
  het_str <- format(het, nsmall = 1)

  prep <- read.csv(
    paste0(RESULTS_FOLDER, "estimates_", het_str, ".csv")
  ) %>%
    preprocess()
  agg <- prep %>%
    group_by(type, T, N) %>%
    summarise(
      fails = TRIALS - n()
    )

  faulty <- unique(agg[agg$fails > 0,]$type)
  cut <- agg[agg$type %in% faulty,]
  cut$rate <- cut$fails / TRIALS

  # Plot the failure rates for different population sizes
  ggplot(
    cut,
    mapping = aes(x = T, y = rate, group = type)
  ) +
    geom_point() +
    geom_line() +
    facet_wrap(
      N ~ type,
      scales = "free_y",
      labeller = labeller(N = label_both, type = label_value)
    ) +
    theme_pubr() +
    ylab("Failure rate")

  ggsave(
    paste0(OUTPUT_FOLDER, "fail_rates_", het_str, ".pdf"),
    width = 210,
    height = 210,
    units = "mm"
  )
}