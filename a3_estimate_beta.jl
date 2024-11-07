include("_100_utils.jl")
include("_110_beta_estimator.jl")
include("_130_benchmarks.jl")

import .BetaEstimator
import .Utils
import .Benchmarks

using StatsBase

import Random: seed!

ALPHAS::Vector{Float64} = [0.5, 2.0]
DATA_FOLDER::String = "./_100_input/simulated/"
breaks_T::Vector{Int64} = collect(5:5:50)
OUTPUT_FOLDER::String = "./_900_output/data/appendix/"
pops::Vector{Int64} = [1000, 5000, 10000]
SEED::Int = 777


seed!(SEED)
for i in eachindex(ALPHAS)
    alpha = ALPHAS[i]
    output_file = OUTPUT_FOLDER * "estimates_$(alpha).csv"
    Utils.write_row(output_file, ["a_hat", "Nu_hat", "N_hat", "No", "trial", "T", "alpha", "N", "type"])
    for N in pops
        data_folder = DATA_FOLDER * "alpha_$(alpha)/"
        data_files = [file for file in readdir(data_folder) if occursin("sample", file) && occursin("$(N).csv", file)]
        for file in data_files[1:100]
            trial_no = parse(Int, split(split(file, "_")[2], ".")[1])
            samples = Utils.read_captures(data_folder * file)
            draws = Array{Any}(undef, 10)
            Threads.@threads for i in eachindex(breaks_T)
                t = breaks_T[i]
                S = samples[1:t]
                O = Set([i for j in S for i in j])
                X = Dict(i => [i in s for s in S] for i in O)
                X = Matrix(
                    transpose(hcat(values(X)...))
                )
                println("***TRIAL NO $file, $t***")
                f = Utils.freq_of_freq(Utils.cap_freq(S))
                initial_N = Benchmarks.turing(length(O), f, t)                
                (minf, minx) = BetaEstimator.fit_Beta(
                    [0.0, log(initial_N - size(X)[1])],
                    X;
                    draws = 1000
                )
                N_hat = length(O) + exp(minx[2])
                println([exp(minx[1]), exp(minx[2]), N_hat, length(O), trial_no, t, alpha, N])
                draws[i] = [exp(minx[1]), exp(minx[2]), N_hat, length(O), trial_no, t, alpha, N, "Beta-binomial"]
            end
            for d in draws
                Utils.write_row(output_file, d)
            end
        end
    end
end
