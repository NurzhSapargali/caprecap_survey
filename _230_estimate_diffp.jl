include("_110_estimator.jl")
include("_100_utils.jl")
include("_120_benchmarks.jl")

import .Estimator
import .Utils
import .Benchmarks

using StatsBase
using DelimitedFiles

import Random: seed!

ALPHAS::Vector{Float64} = [0.5, 1.0, 5.0, 10.0]
DATA_FOLDER::String = "./_200_input/diffp/"
breaks_T::Vector{Int64} = [5, 10, 15, 20]
OUTPUT_FOLDER::String = "./_900_output/data/diffp/"
pops::Vector{Int64} = [1000, 5000, 10000]

function get_truth(filename)
    data = readdlm(filename, ' ', String, '\n')[2, 1]
    return parse(Float64, split(data, ",")[1])
end

for alpha in ALPHAS
    output_file = OUTPUT_FOLDER * "estimates_$(alpha).csv"
    Utils.write_row(output_file, ["a_hat", "Nu_hat", "N_hat", "No", "trial", "T", "alpha", "N", "type"])
    for N in pops
        data_folder = DATA_FOLDER * "alpha_$(alpha)/"
        data_files = [file for file in readdir(data_folder) if occursin("sample", file) && occursin("$(N).csv", file)]
        #N = get_truth(data_folder * "metadata_$(alpha).csv")
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
                No = length(O)
                sum_x = collect(values(K))
                theta = [log(No), log(5.0)]
                theta = Estimator.fit_model(theta, n, No, sum_x)
                N_hat = exp(theta[1]) + No
                println([exp(theta[2]), exp(theta[1]), N_hat, length(O), trial_no, t, alpha, N])
                push!(draws, [exp(theta[2]), exp(theta[1]), N_hat, length(O), trial_no, t, alpha, N, "Pseudolikelihood"])
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
                    push!(draws, [-999.0, benchmarks[b] - No, benchmarks[b], length(O), trial_no, t, alpha, N, b])
                end
            end
            for d in draws
                Utils.write_row(output_file, d)
            end
        end
    end
end
