include("_110_estimator.jl")
include("_100_utils.jl")
include("_120_benchmarks.jl")

import .Estimator
import .Utils
import .Benchmarks

using StatsBase
using DelimitedFiles

BA_EDGES_PER_NODE::Vector{Int64} = [1, 2, 3]
TRIALS::Int64 = 500
breaks_T::Vector{Int64} = [5, 10, 15, 20]
DATA_FOLDER::String = "./_200_input/graphs/"
OUTPUT_FOLDER::String = "./_900_output/data/graphs/"
pops::Vector{Int64} = [1000, 5000, 10000]

function get_truth(filename)
    data = readdlm(filename, ' ', String, '\n')[2, 1]
    return parse(Float64, split(data, ",")[1])
end


for N in pops
    for epn in BA_EDGES_PER_NODE
        output_file = OUTPUT_FOLDER * "nodes_estimates_$(epn).csv"
        Utils.write_row(output_file, ["a_hat", "Nu_hat", "N_hat", "No", "trial", "T", "avg_n", "N", "type"])
        metafile = DATA_FOLDER * "ba_$(epn)/metadata_$N.csv"
        # N = get_truth(metafile)
        data_folder = DATA_FOLDER * "ba_$(epn)/"
        data_files = [file for file in readdir(data_folder) if occursin("nodes", file)]
        for file in data_files
            samples = Utils.read_captures(data_folder * file)
            trial_no = parse(Int, split(split(file, "_")[2], ".")[1])
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
                No = length(O)
                sum_x = collect(values(K))
                theta = [log(No), log(5.0)]
                theta = Estimator.fit_model(theta, n, No, sum_x)
                N_hat = exp(theta[1]) + No
                println([exp(theta[2]), exp(theta[1]), N_hat, length(O), trial_no, t, mean(n), N])
                push!(draws, [exp(theta[2]), exp(theta[1]), N_hat, length(O), trial_no, t, mean(n), N, "Pseudolikelihood"])
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
                    push!(draws, [-999.0, benchmarks[b] - No, benchmarks[b], length(O), trial_no, t, mean(n), N, b])
                end
            end
            for d in draws
                Utils.write_row(output_file, d)
            end
        end
    end
end
