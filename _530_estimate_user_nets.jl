include("./_100_utils.jl")
include("./_120_gamma_estimator.jl")
include("./_130_benchmarks.jl")

import .GammaEstimator
import .Benchmarks
import .Utils

using DataFrames
using CSV

input = "./_900_output/data/user_nets/"

function read_user_net(file)
    df = DataFrame(CSV.File(file))
    return Set(vcat(df[:, :i], df[:, :j]))
end

S = []
for number in 980:986
    file = input * "user_net_$number.csv"
    capture = read_user_net(file)
    push!(S, capture)
end