#= _240_estimate_graphs.jl
This script is designed to estimate population sizes for different sampled
graphs using various methods. It includes utilities for reading data,
benchmarking, and fitting. The script is designed to be run in parallel and
saves the results to CSV files.

Author: Nurzhan Sapargali
Date created: 2024-02-24
=#

# Include necessary files
include("_100_utils.jl")
include("_120_gamma_estimator.jl")
include("_130_benchmarks.jl")

# Import modules
import .GammaEstimator
import .Utils
import .Benchmarks
import Random: seed!

using StatsBase

# Define constants
const GRAPHS::Vector{String} = ["ba_0.05", "ba_0.1", "ba_0.2"] # Folder names containing the sampled graphs
const DATA_FOLDER::String = "./_200_input/graphs/" # Folder containing the sampled graphs
const breaks_T::Vector{Int64} = collect(5:5:50) # Number of samples at which to estimate
const OUTPUT_FOLDER::String = "./_900_output/data/graphs/" # Folder to save the results
const SEED::Int =  777
const TRUE_N::Int =  10000


# Set seed
seed!(SEED)

for i in eachindex(GRAPHS)
    graph = GRAPHS[i]
    output_file = OUTPUT_FOLDER * "estimates_$(graph).csv"

    Utils.write_row(
        output_file,
        ["a_hat", "Nu_hat", "N_hat", "No", "trial", "T", "graph", "N", "type"]
    )

    data_folder = DATA_FOLDER * "$(graph)/"
    data_files = [file for file in readdir(data_folder) if occursin("graph", file)]

    for file in data_files
        trial_no = parse(Int, split(split(file, "_")[2], ".")[1])
        samples = Utils.read_captures(data_folder * file)
        draws = fill([], length(breaks_T) *  15)

        Threads.@threads for j in eachindex(breaks_T)
            t = breaks_T[j]
            S = samples[1:t]
            K = Utils.cap_freq(S)
            f = Utils.freq_of_freq(K)

            println("***TRIAL NO $file, $t***")

            O = Set([i for j in S for i in j])
            n = [length(s) for s in S]
            No = length(O)

            benchmarks = Dict{}()
            benchmarks["Turing"] = Benchmarks.turing(No, f, t)
            benchmarks["Chao"] = Benchmarks.chao(No, f)

            (minf, minx, ret) = GammaEstimator.fit_Gamma([5.0, No], n, No, K)

            N_hat = No + minx[2]

            println([minx[1], minx[2], N_hat, No, trial_no, t, graph,  10000])

            counter =  1
            draws[(j -  1) *  15 + counter] = [
                minx[1],
                minx[2],
                N_hat,
                No,
                trial_no,
                t,
                graph,
                TRUE_N,
                "Pseudolikelihood"
            ]

            # Additional benchmarks
            benchmarks["Schnabel"] = Benchmarks.schnabel(S, n)
            benchmarks["Zelterman"] = Benchmarks.zelterman(No, f)
            benchmarks["Conway-Maxwell-Poisson"] = Benchmarks.conway_maxwell(No, f)
            benchmarks["Turing Geometric"] = Benchmarks.turing_geometric(No, f, t)

            for b in 0:2
                benchmarks["Chao Lee Jeng $b"] = Benchmarks.chao_lee_jeng(
                    No,
                    f,
                    t,
                    n,
                    b
                )
            end

            for k in  1:5
                jk = Benchmarks.jackknife(No, t, f, k)
                benchmarks["Jackknife k = $(k)"] = jk
            end

            for b in keys(benchmarks)
                counter +=  1
                draws[(j -  1) *  15 + counter] = [
                    -999.0,
                    benchmarks[b] - No,
                    benchmarks[b], No,
                    trial_no,
                    t,
                    graph,
                    TRUE_N,
                    b
                ]
            end
        end

        for t in breaks_T
            S = samples[1:t]
            K = Utils.cap_freq(S)
            f = Utils.freq_of_freq(K)
            O = Set([i for j in S for i in j])
            No = length(O)

            mr_estimator = Benchmarks.morgan_ridout(f, t, "./estimateN.R")

            row = [
                -999.0,
                mr_estimator - No,
                mr_estimator,
                No,
                trial_no,
                t,
                graph,
                TRUE_N,
                "Morgan Ridout"
            ]
            push!(draws, row)
        end

        for d in draws
            Utils.write_row(output_file, d)
        end
    end
end
