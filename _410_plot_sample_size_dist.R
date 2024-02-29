library(tidyverse)
library(directlabels)
library(ggpubr)

DENSITIES <- c(0.025, 0.05, 0.1)
INPUT_FOLDER <- "_200_input/graphs/"
OUTPUT_FOLDER <- "./_900_output/figures/" # Folder to save the figures
TRIALS <- 1000

count_sample_sizes <- function(file) {
    res <- readLines(file)
    sapply(gregexpr(",", res), function(x) { length(x) + 1 })
}

label_func <- function(variable_value) {
    return(paste0("Graph density = ", variable_value))
}

sample_sizes <- list()

for (dense in DENSITIES) {
    folder <- paste0(INPUT_FOLDER, "ba_", dense, "/")
    sample_sizes[[as.character(dense)]] <- c()
    for (trial in 1:TRIALS) {
        file <- paste0(folder, "graph_", trial, ".csv")
        sample_sizes[[as.character(dense)]] <- c(
            sample_sizes[[as.character(dense)]],
            count_sample_sizes(file)
        )
    }
}

sizes_df <- do.call(cbind.data.frame, sample_sizes)
cols <- colnames(sizes_df)
sizes_df <- sizes_df %>%
    pivot_longer(
        cols = cols,
        names_to = "density",
        values_to = "sample_size"
    )
ggplot(sizes_df, mapping = aes(x = sample_size)) +
    geom_histogram(fill = "gray", colour = "black", binwidth = 8) +
    theme_minimal() +
    theme_pubr() +
    xlab("Sample size") +
    ylab("Frequency") +
    facet_wrap(
        vars(density),
        ncol = 3,
        labeller = labeller(density = as_labeller(label_func))
    )
ggsave(
    paste0(OUTPUT_FOLDER, "graphs/sample_sizes.pdf"),
    width = 297,
    height = 165,
    units = "mm"
)
