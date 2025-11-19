"""
Simulate capture-recapture data under unequal probability sampling
using a Negative Binomial distribution for sample sizes.
"""

include("_140_simulation_functions.jl")

import .SimFunctions

import Random: seed!

POPS::Vector{Int64} = [1000, 5000, 10000] # Population sizes
DRAWS::Int = 50 # Number of capture occasions
T_MIN::Int = 5 # Minimum number of capture occasions with non-zero recaptures
TRIALS::Int = 1000 # Number of trials per setting
ALPHAS::Vector{Float64} = [0.5, 2.0] # Heterogeneity parameters
NEG_BIN_PARAMS::Vector{Float64} = [1.0, 0.03] # (r, q) pairs for Negative Binomial distribution of sample sizes
SEED::Int = 777 # Random seed for reproducibility
DATA_FOLDER::String = "./_100_input/simulated/" # Folder to save simulated data

# Set random seed
seed!(SEED)

SimFunctions.simulate_data(
    POPS,
    DRAWS,
    ALPHAS,
    NEG_BIN_PARAMS[1],
    NEG_BIN_PARAMS[2],
    DATA_FOLDER,
    TRIALS;
    min_draws = T_MIN
)