include("_100_utils.jl")

using .Utils

import Distributions: Beta, Poisson

import Random: seed!

N::Int = 1000
T::Int = 20
T_min::Int = 5
TRIALS::Int = 5000
ALPHAS::Vector{Float64} = [0.5, 1.0, 5.0, 10.0]
AVG_SAMPLE_SIZE::Int = 30
SEED::Int = 10000
DATA_FOLDER::String = "./_200_input/diffp/";


seed!(SEED)
for alpha in ALPHAS
    metafile = DATA_FOLDER * "alpha_$(alpha)/metadata_$(alpha).csv"
    write_row(metafile, vcat(["N", "T", "alpha"], ["n_$i" for i in 1:T]))
    d = Beta(alpha, (N - 1.0) * alpha)
    p = rand(d, N)
    n = rand(Poisson(AVG_SAMPLE_SIZE), T)
    write_row(metafile, vcat([N, T, alpha], n))
    write_row(metafile, vcat(["p"], p))
    for trial in 1:TRIALS
        println("Generating samples for $n")
        samples = repeat([[]], T)
        O = Set{Int64}()
        while !(length(O) < sum([length(s) for s in samples[1:T_min]]))
            Threads.@threads for i = 1:T
                samples[i] = sampford_sample(p, n[i])
            end
            O = Set([i for j in samples[1:T_min] for i in j])
        end
        file = DATA_FOLDER * "alpha_$(alpha)/sample_$(trial).csv"
        for s in samples
            write_row(file, s)
        end
    end
end
