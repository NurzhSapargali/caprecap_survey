"""
Run estimation methods on simulated data and save results to CSV files.
"""

import Pkg
Pkg.build("RCall")

include("_100_utils.jl")
include("_120_one_nbin.jl")
include("_130_benchmarks.jl")

import .OneNbin
import .Utils
import .Benchmarks

using StatsBase
using Optim

import Random: seed!

ALPHAS::Vector{Float64} = [0.5, 2.0] # Heterogeneity parameters
DATA_FOLDER::String = "./_100_input/simulated/" # Folder with simulated data
breaks_T::Vector{Int64} = collect(5:5:50) # Number of capture occasions to use for estimation
OUTPUT_FOLDER::String = "./_900_output/data/simulated/" # Folder to save estimation results
pops::Vector{Int64} = [1000, 5000, 10000] # Population sizes
N_BENCHMARKS::Int = 0#14 # Number of benchmark methods (excluding Morgan): Schnabel, Chao, Zelterman, CMP, Turing, Turing-G, Chao Lee Jeng (3), Jackknife (5)
SEED::Int = 777

seed!(SEED) # Set random seed for reproducibility
total_estimators = N_BENCHMARKS + 1 # +2 for MPLE-NB and MPLE-G

for alpha in ALPHAS
    # Create output folder and write header to output file
    output_file = OUTPUT_FOLDER * "estimates_$(alpha)_mplenb.csv"
    
    Utils.create_folder_if_not_exists(OUTPUT_FOLDER)
    
    io = open(output_file, "w")
    println(io, "w_hat,a_hat,Nu_hat,N_hat,No,trial,T,alpha,N,type")
    close(io)

    # Read data files for given alpha value 
    data_folder = DATA_FOLDER * "alpha_$(alpha)/"
    for N in pops
        # Filter data files for current population size
        data_files = [
            file for file in readdir(data_folder) if occursin("sample", file) && occursin("$(N).csv", file)
        ]
        for file in data_files
            # Parse trial number from filename and read samples
            trial_no = parse(Int, split(split(file, "_")[2], ".")[1])
            samples = Utils.read_captures(data_folder * file)

            # Initialize array to store results from different methods and numbers of capture occasions
            draws = Array{Any}(nothing, total_estimators * length(breaks_T))

            # Run estimation methods for different numbers of capture occasions
            Threads.@threads for i in eachindex(breaks_T)
                j = 1 # Index for storing results in draws array

                t = breaks_T[i] # Number of capture occasions to use

                S = samples[1:t]
                n = [length(s) for s in S] # Sample sizes

                K = Utils.cap_freq(S)
                f = Utils.freq_of_freq(K)
                No = sum(values(f))

                println("***TRIAL NO $file, $t***")

                # Create dictionary to store benchmark results
                benchmarks = Dict{}()
                benchmarks["Chao"] = Benchmarks.chao(f) # Chao estimator first to use in MPLE initial value

		        initial_N = benchmarks["Chao"] < Inf ? benchmarks["Chao"] : 2 * No
                # Fit MPLE-NB model and store results
                (minf, minx) = OneNbin.fit_oi_nbin_trunc(
                    [log(1.0), log(initial_N - No)],
                    f;
                    upper = [20.0, 23.0],
                    verbose = false,
                    method = Optim.GradientDescent()
                )

                N_hat = No + exp(minx[2])
                w = OneNbin.w_hat(minx[1], minx[2], f)
                a_hat = exp(minx[1])

                for alt_init in [log(No * 3), log(initial_N)]
                    (minf_alt, minx_alt) = OneNbin.fit_oi_nbin_trunc(
                        [log(1.0), alt_init],
                        f;
                        upper = [20.0, 23.0],
                        verbose = false,
                        method = Optim.GradientDescent()
                    )
                    if No + exp(minx_alt[2]) < N_hat
                        N_hat = No + exp(minx_alt[2])
                        w = OneNbin.w_hat(minx_alt[1], minx_alt[2], f)
                        a_hat = exp(minx_alt[1])
                    end
                end

                nbin_row = [w, a_hat, N_hat - No, N_hat, No, trial_no, t, alpha, N, "MPLE-NB"]
                println(nbin_row)

                draws[(i - 1) * total_estimators + j] = nbin_row
                j += 1

                # Fit MPLE-G model and store results
                #(minf, minx) = OneNbin.fit_oi_geom_trunc(
                #    [log(initial_N - No)],
                #    f,
                #    upper = [23.0],
                #    verbose = false
                #)

                #N_hat = No + exp(minx[1])
                #w = OneNbin.w_hat(log(1.0), minx[1], f)

                #geom_row = [w, 1.0, exp(minx[1]), N_hat, No, trial_no, t, alpha, N, "MPLE-G"]
                #println(geom_row)

                #draws[(i - 1) * total_estimators + j] = geom_row
                #j += 1

                # Run benchmark methods and store results except Morgan-Ridout (not parallelizable)
                #benchmarks["Turing"] = Benchmarks.turing(f, t)
                #benchmarks["Schnabel"] = Benchmarks.schnabel(S, n)
                #benchmarks["Zelterman"] = Benchmarks.zelterman(f)
                #benchmarks["Conway-Maxwell-Poisson"] = #Benchmarks.conway_maxwell(f)
                #benchmarks["Turing Geometric"] = Benchmarks.turing_geometric(f)
                #for b in 0:2
                #    benchmarks["Chao Lee Jeng $b"] = Benchmarks.chao_lee_jeng(
                #        f, t, n, b
                #    )
                #end
                #for k in 1:5
                #    jk = Benchmarks.jackknife(t, f, k)
                #    benchmarks["Jackknife k = $(k)"] = jk
                #end
                #for b in keys(benchmarks)
                #    draws[(i - 1) * total_estimators + j] = [
                #        -999.0, -999.0, benchmarks[b] - No, benchmarks[b], No, trial_no, t, alpha, N, b
                    #]
                    #j += 1
                #end
            end

            # Run Morgan-Ridout method separately and store results
            #for t in breaks_T
            #    S = samples[1:t]
            #    K = Utils.cap_freq(S)
            #    f = Utils.freq_of_freq(K)
            #    n = [length(s) for s in S]
            #    No = sum(values(f))
            #    mr_hat = Benchmarks.morgan_ridout(f, t, "./estimateN.R")
            #    push!(
            #        draws,
            #        [-999.0, -999.0, mr_hat - No, mr_hat, No, trial_no, t, alpha, N, "Morgan-Ridout"]
            #    )
            #end

            # Write results to output file
            for d in draws
                Utils.write_row(output_file, d)
            end
        end
    end
end
