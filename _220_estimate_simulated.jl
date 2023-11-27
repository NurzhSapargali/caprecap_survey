include("_100_utils.jl")
include("_120_gamma_estimator.jl")
include("_130_benchmarks.jl")

import .GammaEstimator
import .Utils
import .Benchmarks

using StatsBase

import Random: seed!

ALPHAS::Vector{Float64} = [0.5, 1.0, 5.0, 10.0, Inf]
DATA_FOLDER::String = "./_200_input/caprecap_data/"
breaks_T::Vector{Int64} = collect(5:5:50)
OUTPUT_FOLDER::String = "./_900_output/data/simulated/"
pops::Vector{Int64} = [1000, 5000, 10000]
SEED::Int = 777

seed!(SEED)
for i in 1:length(ALPHAS)
    alpha = ALPHAS[i]
    output_file = OUTPUT_FOLDER * "estimates_$(alpha).csv"
    Utils.write_row(output_file, ["a_hat", "Nu_hat", "N_hat", "No", "trial", "T", "alpha", "N", "type"])
    for N in pops
        data_folder = DATA_FOLDER * "alpha_$(alpha)/"
        data_files = [file for file in readdir(data_folder) if occursin("sample", file) && occursin("$(N).csv", file)]
        for file in data_files
            trial_no = parse(Int, split(split(file, "_")[2], ".")[1])
            samples = Utils.read_captures(data_folder * file)
            draws = []
            for t in breaks_T
                S = samples[1:t]
                K = Utils.cap_freq(S)
                f = Utils.freq_of_freq(K)
                println("***TRIAL NO $file, $t***")
                O = Set([i for j in S for i in j])
                n = [length(s) for s in S]
                No = length(O)
                benchmarks = Dict{}()
                benchmarks["Turing"] = Benchmarks.turing(No, f, t)
                (minf, minx, ret) = GammaEstimator.fit_Gamma(
                    [5.0, Benchmarks.turing(No, f, t) - No],
                    n,
                    No,
                    K
                )
                N_hat = No + minx[2]
                println([minx[1], minx[2], N_hat, No, trial_no, t, alpha, N])
                push!(draws, [minx[1], minx[2], N_hat, No, trial_no, t, alpha, N, "Pseudolikelihood"])
                benchmarks["Schnabel"] = Benchmarks.schnabel(S, n)
                benchmarks["Chao"] = Benchmarks.chao(No, f)
                benchmarks["Zelterman"] = Benchmarks.zelterman(No, f)
                benchmarks["Conway-Maxwell-Poisson"] = Benchmarks.conway_maxwell(No, f)
                benchmarks["Turing Geometric"] = Benchmarks.turing_geometric(No, f, t)
                benchmarks["Morgan Ridout"] = Benchmarks.morgan_ridout(f, t, "./estimateN.R")
                for b in 0:2
                    benchmarks["Chao Lee Jeng $b"] = Benchmarks.chao_lee_jeng(No, f, t, n, b)
                end
                for k in 1:5
                    jk = Benchmarks.jackknife(No, t, f, k)
                    benchmarks["Jackknife k = $(k)"] = jk
                end
                for b in keys(benchmarks)
                    push!(draws, [-999.0, benchmarks[b] - No, benchmarks[b], No, trial_no, t, alpha, N, b])
                end
            end
            for d in draws
                Utils.write_row(output_file, d)
            end
        end
    end
end
