"""
Estimate population size from simulated capture-recapture data under low sample sizes
"""

include("_140_simulation_functions.jl")

import .SimFunctions

import Random: seed!

using StatsBase

ALPHAS::Vector{Float64} = [0.5, 2.0] # Heterogeneity parameters
DATA_FOLDER::String = "./_100_input/simulated/" # Folder with simulated data
BREAKS_T::Vector{Int64} = collect(5:5:50) # Number of capture occasions to use for estimation
OUTPUT_FOLDER::String = "./_900_output/data/appendix/low_sample/" # Folder to save estimation results
POPS::Vector{Int64} = [600, 1000, 5000] # Population sizes
SEED::Int = 777
TRIALS::Int = 1000
INTERMEDIATE_COUNT::Int = 50

 # Indices of files to process for intermediate results to compare against full simulation results
 # Note that the random seed is set below, so each time different files are selected
 # Possible to set specific indices or move this line below the seed!() call for reproducibility
intermediate = sample(1:TRIALS, INTERMEDIATE_COUNT, replace = false)

seed!(SEED) # Set random seed for reproducibility

SimFunctions.estimate_simulations(
    DATA_FOLDER,
    OUTPUT_FOLDER,
    POPS,
    BREAKS_T,
    ALPHAS;
    intermediate = intermediate, # Setting to [0] would run full simulation
    subfolder_suffix = "_low_sample"
)
