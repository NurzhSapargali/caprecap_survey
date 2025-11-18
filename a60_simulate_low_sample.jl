include("_140_simulation_functions.jl")

import .SimFunctions

import Random: seed!

using CSV
using DataFrames

POPS::Vector{Int64} = [600, 1000, 5000] # Population sizes
DRAWS::Int = 50 # Number of capture occasions
T_MIN::Int = 5 # Minimum number of capture occasions with non-zero recaptures
TRIALS::Int = 1000 # Number of trials per setting
ALPHAS::Vector{Float64} = [0.5, 2.0] # Heterogeneity parameters
MEAN_SAMPLE_SIZE::Float64 = 20.0 # Average sample size
SEED::Int = 777 # Random seed for reproducibility
DATA_FOLDER::String = "./_100_input/simulated/" # Folder to save simulated data

# Set random seed
seed!(SEED)

# Calculate the coefficient of variation from existing sampling effort data
sample_efforts = DataFrame(CSV.File("./_100_input/datasets/sampling_effort.csv"))
coef_vars = sample_efforts[:, :coef_var_n]
filter!(!ismissing, coef_vars)
mean_cv = round(sum(coef_vars) / length(coef_vars), digits = 3)

# Parameters of the negative binomial distribution for random sample sizes
var_sample_size = (mean_cv * MEAN_SAMPLE_SIZE)^2.0
q = (MEAN_SAMPLE_SIZE - 1.0) / var_sample_size
r = (MEAN_SAMPLE_SIZE - 1.0)^2.0 / (var_sample_size - MEAN_SAMPLE_SIZE + 1.0)

SimFunctions.simulate_data(
    POPS,
    DRAWS,
    ALPHAS,
    r,
    q,
    DATA_FOLDER,
    TRIALS;
    min_draws = T_MIN,
    subfolder_suffix = "_low_sample"
)
