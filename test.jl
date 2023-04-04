include("_110_estimator.jl")
include("_100_utils.jl")
include("_120_benchmarks.jl")

using .Estimator
using .Utils
using .Benchmarks

using StatsBase
using DelimitedFiles

ALPHAS::Vector{Float64} = [0.5, 1.0, 5.0, 10.0]
DATA_FOLDER::String = "./_200_input/diffp/"
breaks_T::Vector{Int64} = [5, 10, 15, 20]
OUTPUT_FOLDER::String = "./_900_output/data/diffp/"


function get_truth(filename)
    data = readdlm(filename, ' ', String, '\n')[2, 1]
    return parse(Float64, split(data, ",")[1])
end


alpha = 0.5
data_folder = DATA_FOLDER * "alpha_$(alpha)/"
data_files = [file for file in readdir(data_folder) if occursin("sample", file)]
N = get_truth(data_folder * "metadata_$(alpha).csv")
file = "perfect_sample_1.csv"
samples = read_captures(data_folder * file)
t = 5
S = samples[1:t]
K = Dict{Int, Int}()
for s in S
    addcounts!(K, s)
end
f = countmap(values(K))
println("***TRIAL NO $file, $t***")
O = Set([i for j in S for i in j])
n = [length(s) for s in S]
X = Dict{Any, Vector{Bool}}()
for i in O
    X[i] = [i in s for s in S]
    println("....$(length(O) - length(X)) left")
end
(minf, minx, ret) = fit_model(X, n)
