# Plot simulation results for different population size estimators
# using predefined plotting functions under low sample size scenario
source("./_150_plot_functions.R")

RESULTS_FOLDER <- "./_900_output/data/appendix/low_sample/" # Folder containing the simulation results
OUTPUT_FOLDER <- "./_900_output/figures/appendix/" # Folder to save the figures
TRIALS <- 1000 # Number of trials in the simulation
ALPHAS <- c(0.5, 2.0) # Different values of alpha in the simulation
POP_SIZES <- c(600, 1000, 5000) # Different population sizes in the simulation
INTERMEDIATE <- TRUE # Whether to consider file with intermediate results or final results

filename_suffix <- "_low_sample"
if (INTERMEDIATE){
  filename_suffix <- "_low_sample_intermediate"
}

# Plots figures 12, 13, 14, 15 and 16 of the supplementary material corresponding to files:
# ------------------------------
# estimates_box_2.0_low_sample.pdf = Figure 12
# estimates_box_0.5_low_sample.pdf = Figure 13
# relative_bias_low_sample.pdf = Figure 14
# estimates_2.0_low_sample.pdf = Figure 15
# estimates_0.5_low_sample.pdf = Figure 16
# ------------------------------
# in the OUTPUT_FOLDER
plot_everything(
  ALPHAS,
  POP_SIZES,
  NEW_NAMES,
  OLD_NAMES,
  estimates_folder = RESULTS_FOLDER,
  figures_folder = OUTPUT_FOLDER,
  ylim_bias = c(NA, 10),
  ylim_rmse = c(NA, 4.0),
  filename_suffix = filename_suffix, # All output files will have this suffix, e.g., estimates_box_2.0_low_sample_intermediate.pdf
)
