include("_100_utils.jl")

using .Utils

using StatsBase

import Distributions: Beta, Poisson, truncated

import Random: seed!

N::Int = 1000
T::Int = 20
T_min::Int = 5
TRIALS::Int = 500
ALPHAS::Vector{Float64} = [0.5, 1.0, 5.0, 10.0]
n::Vector{Int64} = repeat([5, 100], Int(T / 2))
SEED::Int = 10000
DATA_FOLDER::String = "./_200_input/diffp/"

seed!(SEED)
for alpha in ALPHAS
    metafile = DATA_FOLDER * "alpha_$(alpha)/metadata_$(alpha).csv"
    write_row(metafile, vcat(["N", "T", "alpha"], ["n_$i" for i in 1:T]))
    b = alpha * (N - 1.0)
    d = truncated(Beta(alpha, b), upper = 1.0 / maximum(n))
    p = rand(d, N)
    write_row(metafile, vcat([N, T, alpha], n))
    write_row(metafile, vcat(["p"], p))
    for trial in 1:TRIALS
        println("Generating samples for $n")
        samples = [sample(1:N, pweights(n[i] * p), n[i], replace = false) for i in 1:T]
        file = DATA_FOLDER * "alpha_$(alpha)/sample_$(trial).csv"
        for s in samples
            write_row(file, s)
        end
    end
end
