include("./_100_utils.jl")
include("./_120_gamma_estimator.jl")
include("./_130_benchmarks.jl")

import .GammaEstimator
import .Benchmarks
import .Utils

using StatsBase
using DataFrames
using CSV

const DATA_FOLDER = "./_900_output/data/user_nets/"
const INDICES = [250, 251, 252, 253, 254, 980, 981, 982, 985, 986]

const SAMPLE_SIZES = [233747, 226059, 212918, 192053, 169079]
const PSEUDO = 1.03429e7

function read_user_net(file)
    df = DataFrame(CSV.File(file))
    df = df[df[:, :i] .!= df[:, :j], :]
    return df[:,[:i, :j]]
end

nov_2020 = []
for ind in INDICES
    file = DATA_FOLDER * "user_net_$ind.csv"
    capture = read_user_net(file)
    if ind < 900
        push!(nov_2020, capture)
    end
end

nov_2020 = vcat(nov_2020...)
unique!(nov_2020)
S = []
for i in 1:5
    draw = Set()
    row_inds = collect(1:size(nov_2020, 1))
    n = SAMPLE_SIZES[i] / PSEUDO * size(nov_2020, 1)
    n = round(Int, n)
    remainder = n - length(draw)
    while remainder > 1
        n_edges = Int(round(remainder / 2))
        idx = sample(row_inds, n_edges, replace = false)
        row_inds = setdiff(row_inds, idx)
        union!(draw, nov_2020[idx, :i])
        union!(draw, nov_2020[idx, :j])
        remainder = n - length(draw)
        println(remainder)
    end
    push!(S, draw)
end