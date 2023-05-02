include("_110_estimator.jl")
include("_100_utils.jl")
include("_120_benchmarks.jl")

import .Estimator
import .Utils
import .Benchmarks

using StatsBase
using DelimitedFiles

import Random: seed!

ALPHAS::Vector{Float64} = [1.0]
DATA_FOLDER::String = "./_200_input/diffp/"
breaks_T::Vector{Int64} = [5, 10, 15, 20]
OUTPUT_FOLDER::String = "./_900_output/data/diffp/"
NGRID::Int64 = 100000

function get_truth(filename)
    data = readdlm(filename, ' ', String, '\n')[2, 1]
    return parse(Float64, split(data, ",")[1])
end

for alpha in ALPHAS
    output_file = OUTPUT_FOLDER * "estimates_$(alpha).csv"
    Utils.write_row(output_file, ["a_hat", "b_hat", "N_hat", "No", "trial", "T", "alpha", "N", "type"])
    data_folder = DATA_FOLDER * "alpha_$(alpha)/"
    data_files = [file for file in readdir(data_folder) if occursin("sample", file)]
    N = get_truth(data_folder * "metadata_$(alpha).csv")
    for file in data_files
        trial_no = parse(Int, split(split(file, "_")[2], ".")[1])
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
            (minf, minx, ret) = Estimator.fit_model(X, n, NGRID, [3.0, 6.0 * length(X)]; ftol = 1e-6, upper = [1000.0, Inf])
            N_hat = Estimator.u_size(minx[1], minx[2], maximum(n), length(X)) + length(X)
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
