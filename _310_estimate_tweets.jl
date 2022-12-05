include("_110_estimator.jl")
include("_100_utils.jl")
include("_120_benchmarks.jl")

using .Estimator
using .Utils
using .Benchmarks

using StatsBase
using CSV
using DataFrames

import Random: seed!

MC_DRAWS::Int64 = 1000
DATA_FOLDER::String = "./_900_output/data/hydrated/"

post = [file for file in readdir(DATA_FOLDER) if occursin("_9", file)]
pre = [file for file in readdir(DATA_FOLDER) if occursin("_2", file)]
samples = []
for file in post
    df = DataFrame(CSV.File(DATA_FOLDER * file, delim="\t"))
    net = dropmissing(df, :mentions)
    push!(samples, unique(net[:,:author_id]))
end
K = Dict{Int, Int}()
for s in samples
    addcounts!(K, s)
end
f = countmap(values(K))
O = Set([i for j in samples for i in j])
n = [length(s) for s in samples]
