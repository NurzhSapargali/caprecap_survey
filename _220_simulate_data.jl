include("_100_utils.jl")

using .Utils

using StatsBase

import Distributions: NegativeBinomial, Beta, truncated

import Random: seed!

N::Vector{Int64} = [1000, 5000, 10000]
T::Int = 50
T_min::Int = 5
TRIALS::Int = 1000
ALPHAS::Vector{Float64} = [0.5, 1.0, 5.0, 10.0, Inf]
Q::Float64 = 0.03
R::Int = 1
SEED::Int = 777
DATA_FOLDER::String = "./_200_input/caprecap_data/"

seed!(SEED)
for pop in N
    for alpha in ALPHAS
        metafile = DATA_FOLDER * "alpha_$(alpha)/metadata_$(pop).csv"
        write_row(metafile, ["N", "T", "alpha"])
        write_row(metafile, [pop, T, alpha])
        for trial in 1:TRIALS
            n = rand(NegativeBinomial(R, Q), T) .+ 1
            b = alpha * (pop - 1.0)
            d = truncated(Beta(alpha, b), upper = 1.0 / maximum(n))
            p = rand(d, pop)
            println("Generating samples for $n")
            sum_n_T_min = 0
            samples = []
            O = Set{Int64}()
            if (alpha == Inf)
                # Ensure that there is non-zero recaptures in the first 5 samples
                while !(length(O) < sum_n_T_min)
                    samples = [sample(1:pop, i, replace=false) for i in n]
                    O = Set([i for j in samples[1:T_min] for i in j])
                    sum_n_T_min = sum([length(s) for s in samples[1:T_min]])
                end
            else
                # Ensure that there is non-zero recaptures in the first 5 samples
                while !(length(O) < sum_n_T_min)
                    samples = [sample(1:pop, pweights(i * p), i, replace=false) for i in n]
                    #samples = [sampford_sample(p, i) for i in n]
                    O = Set([i for j in samples[1:T_min] for i in j])
                    sum_n_T_min = sum([length(s) for s in samples[1:T_min]])
                end
            end
            file = DATA_FOLDER * "alpha_$(alpha)/sample_$(trial)_$(pop).csv"
            for s in samples
                write_row(file, s)
            end
        end
    end
end
