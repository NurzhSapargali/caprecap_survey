include("_100_utils.jl")

using .Utils

using StatsBase

import Distributions: Beta, Poisson

import Random: seed!

N::Int = 1000
T::Int = 20
T_min::Int = 5
TRIALS::Int = 1000
ALPHAS::Vector{Float64} = [0.5]
AVG_SAMPLE_SIZE::Int = 30
SEED::Int = 10000
DATA_FOLDER::String = "./_200_input/diffp/"

seed!(SEED)
for alpha in ALPHAS
    metafile = DATA_FOLDER * "alpha_$(alpha)/metadata_perfect_$(alpha).csv"
    write_row(metafile, vcat(["N", "T", "alpha"], ["n_$i" for i in 1:T]))
    d = Beta(alpha, (N - 1.0) * alpha)
    p = rand(d, N)
    n = rand(Poisson(AVG_SAMPLE_SIZE), T)
    write_row(metafile, vcat([N, T, alpha], n))
    write_row(metafile, vcat(["p"], p))
    for trial in 1:TRIALS
        println("Generating samples for $n")
        samples = [sample(1:N, pweights(n[i] * p), n[i], replace = true) for i in 1:T]
        file = DATA_FOLDER * "alpha_$(alpha)/perfect_sample_$(trial).csv"
        for s in samples
            write_row(file, s)
        end
    end
end
