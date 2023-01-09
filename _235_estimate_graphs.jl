include("_110_estimator.jl")
include("_100_utils.jl")
include("_120_benchmarks.jl")

using .Estimator
using .Utils
using .Benchmarks

using StatsBase
using DelimitedFiles

import Random: seed!

BA_EDGES_PER_NODE::Vector{Int64} = [1, 2, 3]
ER_EDGES::Vector{Int64} = [1000, 2000, 3000]
SBM_TYPES::Vector{String} = ["assortative", "disassortative", "core_periphery"]
STRUCTURES::Vector{String} = ["4cliques", "tris", "edges", "nodes"]
TRIALS::Int64 = 100
DATA_FOLDER::String = "./_200_input/graphs/"
OUTPUT_FOLDER::String = "./_900_output/data/graphs/"
MC_DRAWS::Int64 = 2000
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

function read_data(data, unit)
    samples = read_indices(data)
    if (length(samples) <= 1)
        return []
    end
    if (unit == "edges")
        samples = split_indexarray(samples, 2)
    elseif (unit == "tris")
        samples = split_indexarray(samples, 3)
    elseif (unit == "4cliques")
        samples = split_indexarray(samples, 4)
    end
    O = Set([i for j in samples for i in j])
    if sum([length(s) for s in samples]) == length(O)
        return []
    end
    println(samples)
    return samples
end

function get_truth(filename, unit, trial)
    data = readdlm(filename, ' ', Int, '\n')
    if (unit == "edges")
        return data[trial, 2]
    elseif (unit == "tris")
        return data[trial, 3]
    elseif (unit == "4cliques")
        return data[trial, 4]
    end
    return data[trial, 1]
end

function estimate_all(samples, draws, output_dir, trial, truth)
    n = [length(s) for s in samples]
    K = Dict{Any, Int}()
    for s in samples
        for i in s
            if get(K, i, 0) == 0
                K[i] = 1
            else
                K[i] += 1
            end
        end
    end
    f = Dict{Any, Int}()
    for k in values(K)
        if get(f, k, 0) == 0
            f[k] = 1
        else
            f[k] += 1
        end
    end
    t = length(n)
    O = Set([i for j in samples for i in j])
    (minf, minx, ret) = fit_model(samples, O, n, draws)
    v = [minx[1], minx[2] + length(O), minx[2], length(O), trial, length(n), sum(n) / length(n), truth, "Pseudolikelihood"]
    write_row(output_dir, v)
    
    schnab = schnabel(samples, n)
    v[1] = -999
    v[lastindex(v)] = "Schnabel"
    v[2] = schnab
    write_row(output_dir, v)
    
    chao_est = chao(length(O), f)
    v[2] = chao_est
    v[lastindex(v)] = "Chao"
    write_row(output_dir, v)
    
    zelt = zelterman(length(O), f)
    v[2] = zelt
    v[lastindex(v)] = "Zelterman"
    write_row(output_dir, v)
    
    cm = conway_maxwell(length(O), f)
    v[2] = cm
    v[lastindex(v)] = "Conway-Maxwell-Poisson"
    write_row(output_dir, v)
    
    hug = huggins(t, K)
    v[2] = hug
    v[lastindex(v)] = "Huggins"
    write_row(output_dir, v)
    
    alan_geo = turing_geometric(length(O), f, t)
    v[2] = alan_geo
    v[lastindex(v)] = "Turing Geometric"
    write_row(output_dir, v)
    
    alan = turing(length(O), f, t)
    v[2] = alan
    v[lastindex(v)] = "Turing"
    write_row(output_dir, v)
    
    for k in 1:5
        jk = jackknife(length(O), t, f, k)
        v[2] = jk
        v[lastindex(v)] = "Jackknife k = $(k)"
        write_row(output_dir, v)
    end
end

seed!(SEED)
for unit in STRUCTURES
    for trial in 1:TRIALS
        for edges in ER_EDGES
            for gtype in SBM_TYPES
                output = OUTPUT_FOLDER * "sbm_$(edges)/$gtype/$(unit)_estimates.csv"
                filename = DATA_FOLDER * "sbm_$(edges)/$gtype/$(unit)_$(trial).csv"
                S = read_data(filename, unit)
                if (length(S) == 0)
                    continue
                end
                metafile = DATA_FOLDER * "sbm_$(edges)/$gtype/metadata.csv"
                ground_truth = get_truth(metafile, unit, trial)
                println(filename)
                estimate_all(S, MC_DRAWS, output, trial, ground_truth)
            end
            output = OUTPUT_FOLDER * "er_$(edges)/$(unit)_estimates.csv"
            filename = DATA_FOLDER * "er_$(edges)/$(unit)_$(trial).csv"
            S = read_data(filename, unit)
            if (length(S) == 0)
                continue
            end
            metafile = DATA_FOLDER * "er_$(edges)/metadata.csv"
            ground_truth = get_truth(metafile, unit, trial)
            println(filename)
            estimate_all(S, MC_DRAWS, output, trial, ground_truth)
        end
        for epn in BA_EDGES_PER_NODE
            output = OUTPUT_FOLDER * "ba_$(epn)/$(unit)_estimates.csv"
            filename = DATA_FOLDER * "ba_$(epn)/$(unit)_$(trial).csv"
            S = read_data(filename, unit)
            if (length(S) == 0)
                continue
            end
            metafile = DATA_FOLDER * "ba_$(epn)/metadata.csv"
            ground_truth = get_truth(metafile, unit, trial)
            println(filename)
            estimate_all(S, MC_DRAWS, output, trial, ground_truth)
        end
    end
end
