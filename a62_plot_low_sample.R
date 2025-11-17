source("./_150_plot_functions.R")

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

# Plot the results for simulated data
plot_everything(
  ALPHAS,
  POP_SIZES,
  NEW_NAMES,
  OLD_NAMES,
  estimates_folder = RESULTS_FOLDER,
  figures_folder = OUTPUT_FOLDER,
  trials = TRIALS,
  filename_suffix = "_low_sample"
)
