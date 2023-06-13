include("_110_estimator.jl")
include("_100_utils.jl")
include("_120_benchmarks.jl")

import .Estimator
import .Utils
import .Benchmarks

using StatsBase
using DelimitedFiles

import Random: seed!
import StatsFuns: logistic, logit
import Distributions: Normal, truncated

ALPHAS::Vector{Float64} = [0.5]
DATA_FOLDER::String = "./_200_input/diffp/"
breaks_T::Vector{Int64} = [5, 10, 15, 20]
OUTPUT_FOLDER::String = "./_900_output/data/diffp/"
DRAWS::Int64 = 6000
SEED::Int64 = 777

function get_truth(filename)
    data = readdlm(filename, ' ', String, '\n')[2, 1]
    return parse(Float64, split(data, ",")[1])
end

seed!(SEED)
for alpha in ALPHAS
    output_file = OUTPUT_FOLDER * "perfect_estimates_$(alpha).csv"
    #Utils.write_row(output_file, ["a_hat", "b_hat", "N_hat", "No", "trial", "T", "alpha", "N", "type"])
    data_folder = DATA_FOLDER * "alpha_$(alpha)/"
    data_files = [file for file in readdir(data_folder) if occursin("perfect_sample", file)]
    N = get_truth(data_folder * "metadata_perfect_$(alpha).csv")
    #file = "perfect_sample_123.csv"
    for file in data_files[5:500]
        trial_no = parse(Int, split(split(file, "_")[3], ".")[1])
        samples = Utils.read_captures(data_folder * file)
        draws = []
        for t in breaks_T
            S = samples[1:t]
            K = Dict{Int, Int}()
            for s in S
                addcounts!(K, s)
            end
            f = countmap(values(K))
            println("***TRIAL NO $file, $t***")
            O = Set([i for j in S for i in j])
            n = [length(s) for s in S]
            X = Dict{Any, Vector{Bool}}()
            for i in O
                X[i] = [i in s for s in S]
                # println("....$(length(O) - length(X)) left")
            end
            empirical_inclusions = collect(values(K)) / t
            empirical_inclusions[empirical_inclusions .>= 1.0] .= 0.99999
            sigma0 = 0.5 * log(var(logit.(empirical_inclusions)))
            mu0 = logit(1 / length(X))
            (minf, minx, ret) = Estimator.fit_model(X, n, DRAWS, [mu0, sigma0]; ftol = 1e-4)
            mc_draws = logistic.(
                rand(truncated(Normal(minx[1], exp(minx[2])),
                    upper = logit(1 / maximum(n))), 100000)
                )
            N_hat = 1.0 / mean(mc_draws)
#            (minf, minx, ret) = fit_model(X, n, DRAWS, [5.0, length(X)]; upper = [500.0, Inf])
            println([minx[1], minx[2], N_hat, length(O), trial_no, t, alpha, N])
            push!(draws, [minx[1], minx[2], N_hat, length(O), trial_no, t, alpha, N, "Pseudolikelihood"])
            # alpha_trace = [loglh(i, minx[2], S, O, n, 1000) for i in ALPHA_TRACE_RANGE];
            # write_row(OUTPUT_FOLDER * "alpha_trace_$(alpha).csv",
            #           vcat(alpha_trace, [minx[1], minx[2], length(O), t, alpha]));
            # Nu_trace = [loglh(minx[1], i, S, O, n, 1000) for i in NU_TRACE_RANGE];
            # write_row(OUTPUT_FOLDER * "Nu_trace_$(alpha).csv",
            #           vcat(Nu_trace, [minx[1], minx[2], length(O), t, alpha]));
            benchmarks = Dict{}()
            benchmarks["Schnabel"] = Benchmarks.schnabel(S, n)
            benchmarks["Chao"] = Benchmarks.chao(length(O), f)
            benchmarks["Zelterman"] = Benchmarks.zelterman(length(O), f)
            benchmarks["Conway-Maxwell-Poisson"] = Benchmarks.conway_maxwell(length(O), f)
            benchmarks["Huggins"] = Benchmarks.huggins(t, K)
            benchmarks["Turing Geometric"] = Benchmarks.turing_geometric(length(O), f, t)
            benchmarks["Turing"] = Benchmarks.turing(length(O), f, t)
            for k in 1:5
                jk = Benchmarks.jackknife(length(O), t, f, k)
                benchmarks["Jackknife k = $(k)"] = jk
            end
            for b in keys(benchmarks)
                push!(draws, [-999.0, -999.0, benchmarks[b], length(O), trial_no, t, alpha, N, b])
            end
        end
        for d in draws
            Utils.write_row(output_file, d)
        end
    end
end 
