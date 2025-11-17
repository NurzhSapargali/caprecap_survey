"""
Estimate population size from simulated capture-recapture data using various methods.
"""

include("_140_simulation_functions.jl")

using .EstimateSimulations

import Random: seed!

ALPHAS::Vector{Float64} = [0.5, 2.0] # Heterogeneity parameters
DATA_FOLDER::String = "./_100_input/simulated/" # Folder with simulated data
BREAKS_T::Vector{Int64} = collect(5:5:50) # Number of capture occasions to use for estimation
OUTPUT_FOLDER::String = "./_900_output/data/simulated/" # Folder to save estimation results
POPS::Vector{Int64} = [1000, 5000, 10000] # Population sizes
SEED::Int = 777

seed!(SEED) # Set random seed for reproducibility

estimate_simulations(
    DATA_FOLDER,
    OUTPUT_FOLDER,
    POPS,
    BREAKS_T,
    ALPHAS;
    intermediate = false
)
