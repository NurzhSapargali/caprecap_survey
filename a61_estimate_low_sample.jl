"""
Estimate population size from simulated capture-recapture data under low sample sizes
"""

include("_140_simulation_functions.jl")

import .SimFunctions

import Random: seed!

ALPHAS::Vector{Float64} = [0.5, 2.0] # Heterogeneity parameters
DATA_FOLDER::String = "./_100_input/simulated/" # Folder with simulated data
BREAKS_T::Vector{Int64} = collect(5:5:50) # Number of capture occasions to use for estimation
OUTPUT_FOLDER::String = "./_900_output/data/appendix/low_sample/" # Folder to save estimation results
POPS::Vector{Int64} = [600, 1000, 5000] # Population sizes
SEED::Int = 777

seed!(SEED) # Set random seed for reproducibility

SimFunctions.estimate_simulations(
    DATA_FOLDER,
    OUTPUT_FOLDER,
    POPS,
    BREAKS_T,
    ALPHAS;
    #intermediate = true,
    subfolder_suffix = "_low_sample"
)
