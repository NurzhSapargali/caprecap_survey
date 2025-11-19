"""
Functions to simulate unequal probability sampling and run estimation methods 
on simulated data.
"""
module SimFunctions

import Pkg

include("_100_utils.jl")
include("_120_one_nbin.jl")
include("_130_benchmarks.jl")

import .Utils
import .Benchmarks

using Optim
using StatsBase

export estimate_simulations, simulate_data


"""
    simulate_data(
        pops::Vector{Int64},
        draws::Int,
        alphas::Vector{Float64},
        r::Float64,
        q::Float64,
        data_folder::String,
        trials::Int;
        min_draws::Int64 = 5,
        subfolder_suffix::String = ""
    )

Simulate capture–recapture data for given population sizes (`pops`), numbers of
capture occasions (`draws`), heterogeneity parameters (`alphas`), and Negative
Binomial sample size parameters (`r`, `q`). Save the simulated data to `data_folder`.
Each setting is repeated `trials` times. The optional parameter `min_draws` sets
the minimum number of captures per individual. An optional `subfolder_suffix` can
be added to the data subfolder names.
"""
function simulate_data(
    pops::Vector{Int64},
    draws::Int,
    alphas::Vector{Float64},
    r::Float64,
    q::Float64,
    data_folder::String,
    trials::Int;
    min_draws::Int64 = 5,
    subfolder_suffix::String = ""
)
    # Create data folder and subfolders for different alpha values
    for alpha in alphas
        Utils.create_folder_if_not_exists(
            data_folder * "alpha_$(alpha)" * subfolder_suffix * "/"
        )
    end

    # Loop over population sizes, heterogeneity parameters, and sample size parameters
    for N in pops
        for alpha in alphas
            # Write to metadata file current setting parameters
            metafile = data_folder * "alpha_$(alpha)" * subfolder_suffix * "/metadata_$(N).csv"
            Utils.write_row(metafile, ["N", "T", "alpha", "r", "q"])
            Utils.write_row(metafile, [N, draws, alpha, r, q])

            # Repeatedly simulate data and save to files
            for trial in 1:trials
                println(
                    "Simulating data for N=$N, alpha=$alpha, r=$r, q=$q, trial=$trial"
                )
                samples = Utils.simulate_samples(N, draws, alpha, r, q, min_draws)
                println(
                    "Unique captured individuals: ",
                    length(Set([i for s in samples for i in s]))
                )
                println("-----")

                # Save sample lists to file
                file = data_folder * "alpha_$(alpha)" * subfolder_suffix * "/sample_$(trial)_$(N).csv"
                for s in samples
                    Utils.write_row(file, s)
                end
            end
        end
    end
end


"""
    estimate_simulations(
        data_folder::String,
        output_folder::String,
        pops::Vector{Int64},
        breaks_T::Vector{Int64},
        alphas::Vector{Float64};
        intermediate = false,
        subfolder_suffix::String = ""
    )

Run estimation methods on simulated data located in `data_folder` for the specified
population sizes (`pops`), numbers of capture occasions (`breaks_T`), and heterogeneity
parameters (`alphas`). Save the results to CSV files in `output_folder`.
If `intermediate` is true, randomly select 100 data files per setting for estimation.
An optional `subfolder_suffix` can be added to the data subfolder names.
"""
function estimate_simulations(
    data_folder::String,
    output_folder::String,
    pops::Vector{Int64},
    breaks_T::Vector{Int64},
    alphas::Vector{Float64};
    intermediate = false,
    subfolder_suffix::String = ""
)
    Pkg.build("RCall")
    total_estimators = 14 + 2 # benchmarks + 2 for MPLE-NB and MPLE-G

    for alpha in alphas
        # Create output folder and write header to output file
        output_file = output_folder * "estimates_$(alpha)" * subfolder_suffix * ".csv"
        # For intermediate results, set different output file name
        if intermediate
            output_file = output_folder * "estimates_$(alpha)" * subfolder_suffix * "_intermediate.csv"
        end
    
        Utils.create_folder_if_not_exists(output_folder)
    
        io = open(output_file, "w")
        println(io, "w_hat,a_hat,Nu_hat,N_hat,No,trial,T,alpha,N,type")
        close(io)

        # Read data files for given alpha value 
        subfolder = data_folder * "alpha_$(alpha)" * subfolder_suffix * "/"
        for N in pops
            # Filter data files for current population size
            data_files = [
                file for file in readdir(subfolder) if occursin("sample", file) && occursin("$(N).csv", file)
            ]
            # If intermediate flag is set, randomly select 100 files
            if intermediate
                random_index = sample(1:length(data_files), 100; replace = false)
                data_files = data_files[random_index]
            end

            for file in data_files
                # Parse trial number from filename and read samples
                trial_no = parse(Int, split(split(file, "_")[2], ".")[1])
                samples = Utils.read_captures(subfolder * file)

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

                    # Initialize at Chao if finite, else at 2*No
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

                    # Try alternative initializations for better optimum
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

                    nbin_row = [
                        w, a_hat, N_hat - No, N_hat, No, trial_no, t, alpha, N, "MPLE-NB"
                    ]
                    println(nbin_row)

                    draws[(i - 1) * total_estimators + j] = nbin_row
                    j += 1

                    # Fit MPLE-G model and store results
                    (minf, minx) = OneNbin.fit_oi_geom_trunc(
                        [log(initial_N - No)],
                        f,
                        upper = [23.0],
                        verbose = false
                    )

                    N_hat = No + exp(minx[1])
                    w = OneNbin.w_hat(log(1.0), minx[1], f)

                    geom_row = [
                        w, 1.0, exp(minx[1]), N_hat, No, trial_no, t, alpha, N, "MPLE-G"
                    ]
                    println(geom_row)

                    draws[(i - 1) * total_estimators + j] = geom_row
                    j += 1

                    # Run benchmark methods and store results except Morgan-Ridout (not parallelizable)
                    benchmarks["Turing"] = Benchmarks.turing(f, t)
                    benchmarks["Schnabel"] = Benchmarks.schnabel(S, n)
                    benchmarks["Zelterman"] = Benchmarks.zelterman(f)
                    benchmarks["Conway-Maxwell-Poisson"] = #Benchmarks.conway_maxwell(f)
                    benchmarks["Turing Geometric"] = Benchmarks.turing_geometric(f)
                    for b in 0:2
                        benchmarks["Chao Lee Jeng $b"] = Benchmarks.chao_lee_jeng(
                            f, t, n, b
                        )
                    end
                    for k in 1:5
                        jk = Benchmarks.jackknife(t, f, k)
                        benchmarks["Jackknife k = $(k)"] = jk
                    end
                    for b in keys(benchmarks)
                        bench_row = [
                            -999.0,
                            -999.0,
                            benchmarks[b] - No,
                            benchmarks[b],
                            No,
                            trial_no,
                            t,
                            alpha,
                            N,
                            b
                        ]
                        draws[(i - 1) * total_estimators + j] = bench_row
                        j += 1
                    end
                end

                # Run Morgan-Ridout method separately and store results
                for t in breaks_T
                    S = samples[1:t]
                    K = Utils.cap_freq(S)
                    f = Utils.freq_of_freq(K)
                    n = [length(s) for s in S]
                    No = sum(values(f))

                    mr_hat = Benchmarks.morgan_ridout(f, t, "./estimateN.R")
                    mr_row = [
                        -999.0,
                        -999.0,
                        mr_hat - No,
                        mr_hat,
                        No,
                        trial_no,
                        t,
                        alpha,
                        N,
                        "Morgan-Ridout"
                    ]
                    push!(draws, mr_row)
                end

                # Write results to output file
                for d in draws
                    Utils.write_row(output_file, d)
                end
            end
        end
    end
end

end
