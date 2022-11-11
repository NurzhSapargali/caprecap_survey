include("_110_estimator.jl")
include("_100_utils.jl")

using .Estimator
using .Utils

using StatsBase

import Random: seed!

ALPHA::Float64 = 5.0
DATA_FOLDER::String = "./_200_input/diffp/"
T::Int64 = 10
OUTPUT_FOLDER::String = "./_900_output/data/diffp/"
MC_DRAWS::Int64 = 2500
SEED::Int64 = 111
TRIALS::Int64 = 10

function leave_one_out(S, index, draws)
    reduced = copy(S)
    deleteat!(reduced, index)
    O = Set([i for j in reduced for i in j])
    n = [length(s) for s in reduced]
    (minf, minx, ret) = fit_model(reduced, O, n, draws)
    out = copy(minx)
    push!(out, length(O))
    return out
end

seed!(SEED);
folder = DATA_FOLDER * "alpha_$(ALPHA)/";
files = [file for file in readdir(folder) if occursin("sample", file)][10:TRIALS + 10]
res = []
for file in files
    samples = read_captures(folder * file)[1:T]
    jack_estimates = [[0.0, 0.0, 0.0] for i in 1:T]
    Threads.@threads for i = 1:T
        jack_estimates[i] = leave_one_out(samples, i, MC_DRAWS)
    end
    Njacks = [i[2] + i[3] for i in jack_estimates]
    sterr = sqrt((T - 1) / T * sum((Njacks .- mean(Njacks)).^2))
    push!(res, sterr)
end
println(res)
open("jackknife_devs2.csv", "w") do f
    for i in res
        println(f, i)
    end
end
