# Plot simulation results for different population size estimators
# using predefined plotting functions.
source("./_150_plot_functions.R")

RESULTS_FOLDER <- "./_900_output/data/simulated/" # Folder containing the simulation results
OUTPUT_FOLDER <- "./_900_output/figures/simulated/" # Folder to save the figures
ALPHAS <- c(0.5, 2.0) # Different values of alpha in the simulation
POP_SIZES <- c(1000, 5000, 10000) # Different population sizes in the simulation
INTERMEDIATE <- FALSE # Whether to consider file with intermediate results or final results

filename_suffix <- ""
if (INTERMEDIATE){
  filename_suffix <- "_intermediate"
}

# Plots figures 1, 2, 3, 4 and 5 of the main text corresponding to files:
# ------------------------------
# estimates_box_2.0.pdf = Figure 1
# estimates_box_0.5.pdf = Figure 2
# relative_bias.pdf = Figure 3
# estimates_2.0.pdf = Figure 4
# estimates_0.5.pdf = Figure 5
# ------------------------------
# in the OUTPUT_FOLDER
plot_everything(
  ALPHAS,
  POP_SIZES,
  NEW_NAMES,
  OLD_NAMES,
  estimates_folder = RESULTS_FOLDER,
  figures_folder = OUTPUT_FOLDER,
  ylim_bias = c(NA, 7.0),
  ylim_rmse = c(NA, 4.0),
  filename_suffix = filename_suffix # all output files will have this suffix, e.g., estimates_box_2.0_intermediate.pdf
)
