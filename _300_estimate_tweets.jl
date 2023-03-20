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
SEED::Int = 123

seed!(SEED)
post = [file for file in readdir(DATA_FOLDER) if (occursin("_9", file)&occursin(".csv", file))]
pre = [file for file in readdir(DATA_FOLDER) if (occursin("_2", file)&occursin(".csv", file))]
samples = []
networks = []
for file in post
    df = DataFrame(CSV.File(DATA_FOLDER * file, delim="\t"))
    net = dropmissing(df, :mentions)
    self_loop = string.(net[:, :author_id]) .== [i[3:lastindex(i) - 2] for i in net[:, :mentions]]
    net = net[.!(self_loop), :]
    adj_list = []
    for i in 1:size(net)[1]
        cut = net[i, :mentions]
        cut = cut[3:lastindex(cut) - 2]
        adjacents = split(cut, "', '")
        for u in adjacents
            push!(adj_list, (net[i, :author_id],  parse(Int, u)))
        end
    end
    push!(samples, unique([i for j in adj_list for i in j]))
    push!(networks, adj_list)
end

K = Dict{Int, Int}()
for s in samples
    addcounts!(K, s)
end
f = countmap(values(K))
O = Set([i for j in samples for i in j])
n = [length(s) for s in samples]
