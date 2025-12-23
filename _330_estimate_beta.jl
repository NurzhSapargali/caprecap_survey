"""
Estimate population size using Beta-binomial model on PPS capture-recapture data
simulated under different heterogeneity settings. Results are saved in output
CSV files.
"""

include("_100_utils.jl")
include("_110_beta_estimator.jl")
include("_130_benchmarks.jl")

import .BetaEstimator
import .Utils
import .Benchmarks

using StatsBase

import Random: seed!

ALPHAS::Vector{Float64} = [0.5, 2.0] # Heterogeneity parameters
DATA_FOLDER::String = "./_100_input/simulated/" # Folder containing simulated datasets
BREAKS_T::Vector{Int64} = collect(5:5:50) # Numbers of capture occasions to consider
OUTPUT_FOLDER::String = "./_900_output/data/appendix/beta_bin/" # Output folder
POPS::Vector{Int64} = [1000, 5000, 10000] # Population sizes to consider
MC_DRAWS::Int = 1000 # Number of Monte Carlo draws for Beta-binomial estimation
SEED::Int = 777
TRIALS::Int = 100
INTERMEDIATE_COUNT::Int = parse(Int, get(ENV, "JULIA_INTERMEDIATE_COUNT", "0")) # Set to 0 to run full simulation

 # Indices of files to process for intermediate results to compare against full simulation results
 # Note that the random seed is set below, so each time different files are selected
 # Possible to set specific indices or move this line below the seed!() call for reproducibility
intermediate = [] # Default to full simulation
if INTERMEDIATE_COUNT > 0
    intermediate = sample(1:TRIALS, INTERMEDIATE_COUNT, replace = false)
end


seed!(SEED)
u = rand(MC_DRAWS) # Precompute uniform draws for Monte Carlo integration

for i in eachindex(ALPHAS)
    alpha = ALPHAS[i]
    # Create output folder and write header to output file
    output_file = OUTPUT_FOLDER * "estimates_$(alpha)_betabin.csv"
    # If intermediate indices are specified, adjust output filename
    if length(intermediate) > 0
        output_file = OUTPUT_FOLDER * "estimates_$(alpha)_betabin_intermediate.csv"
    end

    Utils.create_folder_if_not_exists(OUTPUT_FOLDER)

    io = open(output_file, "w")
    println(io, "a_hat,Nu_hat,N_hat,No,trial,T,alpha,N,type")
    close(io)

    for N in POPS
        data_folder = DATA_FOLDER * "alpha_$(alpha)/"
        # Process only first 100 files to limit runtime
        data_files = [
            file for file in readdir(data_folder)
            if occursin("sample", file) && occursin("$(N).csv", file)
        ][1:TRIALS]

        # If intermediate indices are specified, select only those files
        if length(intermediate) > 0
            data_files = data_files[intermediate]
        end
        rows = Array{Any}(nothing, length(data_files) * length(BREAKS_T))
        for j in eachindex(data_files)
            file = data_files[j]
            # Parse trial number from filename and read samples
            trial_no = parse(Int, split(split(file, "_")[2], ".")[1])
            samples = Utils.read_captures(data_folder * file)

            # Initialize array to store results from different methods and T
            draws = Array{Any}(nothing, length(BREAKS_T))
            
            # Run estimation methods for different numbers of capture occasions
            println("***TRIAL NO $file***")
            Threads.@threads for i in eachindex(BREAKS_T)
                t = BREAKS_T[i] # Number of capture occasions to use
                S = samples[1:t]

                O = Set([i for j in S for i in j]) # Observed individuals
                X = Dict(i => [i in s for s in S] for i in O) # Capture history matrix
                X = Matrix(
                    transpose(hcat(values(X)...))
                )

                f = Utils.freq_of_freq(Utils.cap_freq(S))

                # Determine initial population size using Turing estimator
                initial_N = Benchmarks.turing(f, t) < Inf ? Benchmarks.turing(f, t) : 2 * length(O)

                # Fit Beta-binomial model and store results
                (minf, minx) = BetaEstimator.fit_Beta(
                    [0.0, log(initial_N - size(X)[1])],
                    X,
                    u;
                    verbose = false
                )

                N_hat = length(O) + exp(minx[2]) # Estimated population size
                println(
                    [exp(minx[1]), exp(minx[2]), N_hat, length(O), trial_no, t, alpha, N]
                )

                draws[i] = [
                    exp(minx[1]), # a_hat
                    exp(minx[2]), # Nu_hat
                    N_hat,
                    length(O),
                    trial_no,
                    t,
                    alpha,
                    N,
                    "Beta-binomial"
                ]
            end
            # Store results in rows array
            for i in eachindex(BREAKS_T)
                Utils.write_row(output_file, draws[i])
            end
        end
        # Write results to output file
    end
end
