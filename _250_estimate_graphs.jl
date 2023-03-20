include("_110_estimator.jl")
include("_100_utils.jl")
include("_120_benchmarks.jl")

using .Estimator
using .Utils
using .Benchmarks

using StatsBase
using DelimitedFiles

BA_EDGES_PER_NODE::Vector{Int64} = [1]
# SBM_TYPES::Vector{String} = ["assortative", "disassortative", "core_periphery"]
TRIALS::Int64 = 5000
DATA_FOLDER::String = "./_200_input/graphs/"
OUTPUT_FOLDER::String = "./_900_output/data/graphs/"
GRID_SIZE::Int64 = 10000

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

function read_data(data)
    samples = read_indices(data)
    if (length(samples) <= 1)
        return []
    end
    O = Set([i for j in samples for i in j])
    if sum([length(s) for s in samples]) == length(O)
        return []
    end
    println(samples)
    return samples
end

function get_truth(filename)
    data = readdlm(filename, ' ', String, '\n')[2, 1]
    return parse(Float64, split(data, ",")[1])
end

function estimate_all(samples, ngrid, output_dir, trial, truth)
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
    (minf, minx, ret) = fit_model(samples, O, n, ngrid)
    v = [minx[1], minx[2] + length(O), minx[2], length(O), trial, t, sum(n) / length(n), truth, "Pseudolikelihood"]
    write_row(output_dir, v)
    benchmarks = Dict{}()
    benchmarks["Schnabel"] = schnabel(samples, n)
    benchmarks["Chao"] = chao(length(O), f)
    benchmarks["Zelterman"] = zelterman(length(O), f)
    benchmarks["Conway-Maxwell-Poisson"] = conway_maxwell(length(O), f)
    benchmarks["Huggins"] = huggins(t, K)
    benchmarks["Turing Geometric"] = turing_geometric(length(O), f, t)
    benchmarks["Turing"] = turing(length(O), f, t)
    for k in 1:5
        jk = jackknife(length(O), t, f, k)
        benchmarks["Jackknife k = $(k)"] = jk
    end
    for b in keys(benchmarks)
        write_row(output_dir,
                  [-999.0, benchmarks[b], minx[2], length(O), trial, t, sum(n) / length(n), truth, b])
    end
end


for epn in BA_EDGES_PER_NODE
    output = OUTPUT_FOLDER * "ba_$(epn)/nodes_estimates.csv"
    write_row(output, ["a_hat", "N_hat", "Nu_hat", "No", "trial", "T", "avg_n", "N", "type"])
    metafile = DATA_FOLDER * "ba_$(epn)/metadata.csv"
    ground_truth = get_truth(metafile)
    for trial in 1:TRIALS
        filename = DATA_FOLDER * "ba_$(epn)/nodes_$(trial).csv"
        S = read_captures(filename)
        if (length(S) == 0)
            continue
        end
        println(filename)
        estimate_all(S, GRID_SIZE, output, trial, ground_truth)
    end
end
