include("_110_estimator.jl")
include("_100_utils.jl")
include("_120_benchmarks.jl")

using .Estimator
using .Utils
using .Benchmarks

using StatsBase
using DelimitedFiles

import Random: seed!
import Distributions: Beta, logcdf

ALPHAS::Vector{Float64} = [0.5]
DATA_FOLDER::String = "./_200_input/diffp/"
breaks_T::Vector{Int64} = [5, 10, 15, 20, 30, 35]
OUTPUT_FOLDER::String = "./_900_output/data/diffp/"
DRAWS::Int = 35000


function get_truth(filename)
    data = readdlm(filename, ' ', String, '\n')[2, 1]
    return parse(Float64, split(data, ",")[1])
end

function logN(a, b, upper)
    main_term = (log(a + b) - log(a))
    temper = logcdf(Beta(a, b), upper) - logcdf(Beta(a + 1, b), upper)
    println(exp(main_term + temper))
    return main_term + temper
end

for alpha in ALPHAS
    output_file = OUTPUT_FOLDER * "perfect_estimates_$(alpha).csv"
    #write_row(output_file, ["a_hat", "N_hat", "Nu_hat", "No", "trial", "T", "alpha", "N", "type"])
    data_folder = DATA_FOLDER * "alpha_$(alpha)/"
    data_files = [file for file in readdir(data_folder) if occursin("perfect_sample", file)]
    N = get_truth(data_folder * "metadata_perfect_$(alpha).csv")
    for file in data_files
        trial_no = parse(Int, split(split(file, "_")[3], ".")[1])
        samples = read_captures(data_folder * file)
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
                println("....$(length(O) - length(X)) left")
            end
            (minf, minx, ret) = fit_model(X, n, [5.0, length(X)]; ftol = 0.000000001)
#            (minf, minx, ret) = fit_model(X, n, DRAWS, [5.0, length(X)]; upper = [500.0, Inf])
            write_row(output_file,
                      [minx[1], exp(logN(minx[1], minx[2], 1 / maximum(n))), minx[2], length(O), trial_no, t, alpha, N, "Pseudolikelihood"])
            # alpha_trace = [loglh(i, minx[2], S, O, n, 1000) for i in ALPHA_TRACE_RANGE];
            # write_row(OUTPUT_FOLDER * "alpha_trace_$(alpha).csv",
            #           vcat(alpha_trace, [minx[1], minx[2], length(O), t, alpha]));
            # Nu_trace = [loglh(minx[1], i, S, O, n, 1000) for i in NU_TRACE_RANGE];
            # write_row(OUTPUT_FOLDER * "Nu_trace_$(alpha).csv",
            #           vcat(Nu_trace, [minx[1], minx[2], length(O), t, alpha]));
            benchmarks = Dict{}()
            benchmarks["Schnabel"] = schnabel(S, n)
            benchmarks["Chao"] = chao(length(O), f)
            benchmarks["Zelterman"] = zelterman(length(O), f)
            benchmarks["Conway-Maxwell-Poisson"] = conway_maxwell(length(O), f)
            benchmarks["Huggins"] = huggins(t, K)
            benchmarks["Turing Geometric"] = turing_geometric(length(O), f, t)
            benchmarks["Turing"] = turing(length(O), f, t)
            for k in 1:5
                jk = jackknife(length(O), t, f, k)
                benchmarks["Jackknife k = $(k)"] = jk
            end
            for b in keys(benchmarks)
                write_row(output_file,
                          [-999.0, benchmarks[b], minx[2], length(O), trial_no, t, alpha, N, b])
            end
        end
    end
end 