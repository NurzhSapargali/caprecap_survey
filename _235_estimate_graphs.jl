include("_110_estimator.jl")
include("_100_utils.jl")

using .Estimator
using .Utils

import Random: seed!

EDGES_PER_NODE::Vector{Int64} = [1, 2, 3]
STRUCTURES::Vector{String} = ["4cliques", "tris", "edges", "nodes"]
TRIALS::Int64 = 25
DATA_FOLDER::String = "./_200_input/graphs/"
OUTPUT_FOLDER::String = "./_900_output/data/graphs/"
MC_DRAWS::Int64 = 5000
SEED::Int64 = 111

function read_indices(filename::String)
    S = try 
        [map(x -> parse(Int, x), split(line, ",")) for line in eachline(filename)]
    catch error
        if isa(error, SystemError)
            []
        end
    end
    return S
end

function split_indexarray(inds::Vector{Vector{Int64}},
                          step::Int64)
    out = []
    for s in inds
        if (step == 2)
            row = [sort([s[i], s[i + 1]]) for i in 1:step:length(s)]
        elseif (step == 3)
            row = [sort([s[i], s[i + 1], s[i + 2]]) for i in 1:step:length(s)]
        elseif (step == 4)
            row = [sort([s[i], s[i + 1], s[i + 2], s[i + 3]]) for i in 1:step:length(s)]
        end
        push!(out, row)
    end
    return out
end


seed!(SEED)
for epn in EDGES_PER_NODE
    for unit in STRUCTURES
        output = OUTPUT_FOLDER * "ba_$(epn)/$(unit)_estimates.csv"
        for trial in 1:TRIALS
            filename = DATA_FOLDER * "ba_$(epn)/$(unit)_$(trial).csv"
            S = read_indices(filename)
            if (length(S) == 0)
                continue
            end
            if (unit == "edges")
                S = split_indexarray(S, 2)
            elseif (unit == "tris")
                S = split_indexarray(S, 3)
            elseif (unit == "4cliques")
                S = split_indexarray(S, 4)
            end
            println(S)
            O = Set([i for j in S for i in j])
            n = [length(s) for s in S]
            (minf, minx, ret) = fit_model(S, O, n, MC_DRAWS);
            write_row(output,
                      [minx[1], minx[2], length(O), length(n), sum(n) / length(n)]);
        end
    end
end
