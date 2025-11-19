# Plot simulation results for different population size estimators
# using predefined plotting functions under low population size scenario
source("./_150_plot_functions.R")

RESULTS_FOLDER <- "./_900_output/data/appendix/low_pop/" # Folder containing the simulation results
OUTPUT_FOLDER <- "./_900_output/figures/appendix/" # Folder to save the figures
ALPHAS <- c(0.5, 2.0) # Different values of alpha in the simulation
POP_SIZES <- c(600, 800) # Different population sizes in the simulation

# Plots figures 7, 8, 9, 10 and 11 of the supplementary material corresponding to files:
# ------------------------------
# estimates_box_2.0_low_pop.pdf = Figure 7
# estimates_box_0.5_low_pop.pdf = Figure 8
# relative_bias_low_pop.pdf = Figure 9
# estimates_2.0_low_pop.pdf = Figure 10
# estimates_0.5_low_pop.pdf = Figure 11
# ------------------------------
# in the OUTPUT_FOLDER
plot_everything(
  ALPHAS,
  POP_SIZES,
  NEW_NAMES,
  OLD_NAMES,
  estimates_folder = RESULTS_FOLDER,
  figures_folder = OUTPUT_FOLDER,
  filename_suffix = "_low_pop",
  box_facet_size = c(198, 210),
  comp_facet_size = c(198, 210)
)

# Plots figures from intermediate data (random 100 simulated datasets per setting)
# The filenames are as above but with suffix "_intermediate", e.g.,
# estimates_box_2.0_low_pop_intermediate.pdf, estimates_0.5_low_pop_intermediate.pdf
plot_everything(
  ALPHAS,
  POP_SIZES,
  NEW_NAMES,
  OLD_NAMES,
  estimates_folder = RESULTS_FOLDER,
  figures_folder = OUTPUT_FOLDER,
  filename_suffix = "_low_pop_intermediate",
  ylim_bias = c(NA, 8.0), # Approximately same limits as with all simulated data
  ylim_rmse = c(NA, 4.0),
  box_facet_size = c(198, 210),
  comp_facet_size = c(198, 210)
)
