include("_100_utils.jl")

using .Utils

import Random: seed!

N::Vector{Int64} = [600, 800] # Population sizes
T::Int = 50 # Number of capture occasions
T_min::Int = 5 # Minimum number of capture occasions with non-zero recaptures
TRIALS::Int = 1000 # Number of trials per setting
ALPHAS::Vector{Float64} = [0.5, 2.0] # Heterogeneity parameters
NEG_BIN_PARAMS::Vector{Float64} = [1.0, 0.03] # (r, q) pairs for Negative Binomial distribution of sample sizes
SEED::Int = 777 # Random seed for reproducibility
DATA_FOLDER::String = "./_100_input/simulated/" # Folder to save simulated data

# Set random seed
seed!(SEED)

# Create data folder and subfolders for different alpha values
Utils.create_folder_if_not_exists(DATA_FOLDER)
for alpha in ALPHAS
    Utils.create_folder_if_not_exists(DATA_FOLDER * "alpha_$(alpha)_low_pop/")
end

r = NEG_BIN_PARAMS[1]
q = NEG_BIN_PARAMS[2]
# Loop over population sizes, heterogeneity parameters, and sample size parameters
for pop in N
    for alpha in ALPHAS

        # Write to metadata file current setting parameters
        metafile = DATA_FOLDER * "alpha_$(alpha)_low_pop/metadata_$(pop).csv"
        Utils.write_row(metafile, ["N", "T", "alpha", "r", "q"])
        Utils.write_row(metafile, [pop, T, alpha, r, q])

        # Repeatedly simulate data and save to files
        for trial in 1:TRIALS
            println("Simulating data for N=$pop, alpha=$alpha, r=$r, q=$q, trial=$trial")
            samples = Utils.simulate_samples(pop, T, alpha, r, q)
            println("Unique captured individuals: ", length(Set([i for s in samples for i in s])))
            println("-----")

            # Save sample lists to file
            file = DATA_FOLDER * "alpha_$(alpha)_low_pop/sample_$(trial)_$(pop).csv"
            for s in samples
                Utils.write_row(file, s)
            end
        end
    end
end
