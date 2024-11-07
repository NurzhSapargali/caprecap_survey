include("_100_utils.jl")
include("_120_one_nbin.jl")
include("_130_benchmarks.jl")

import .OneNbin
import .Utils
import .Benchmarks

using StatsBase
using CSV
using DataFrames

import Random: seed!

INPUT_FOLDER::String = "./_100_input/datasets/"
OUTPUT_FOLDER::String = "./_900_output/data/datasets/"
SEED::Int = 777

seed!(SEED)
data_files = [file for file in readdir(INPUT_FOLDER) if file != "metadata.csv"]
df = DataFrame(
    dataset = String[],
    No = Int64[],
    N_hat_nb = Float64[],
    w_hat_nb = Float64[],
    a_hat = Float64[],
    N_hat_geom = Float64[],
    w_hat_geom = Float64[]
)
for file in data_files
    println("***$file***")
    file_path = joinpath(INPUT_FOLDER, file)
    new_data = CSV.read(file_path, DataFrame, header=false)
    f = Dict(new_data[:, 1] .=> new_data[:, 2])
    No = sum(values(f))
    chao = Benchmarks.chao(No, f)
    nb_cands = Dict{}()
    geom_cands = Dict{}()
    for initial in [log(2.0), log(chao - No), log(No)]                   
        (minf, minx) = OneNbin.fit_oi_nbin_trunc(
            [0.0, log(1.0), initial],
            f,
            upper = [Inf, 20.0, 23.0],
            verbose = false
        )
        nb_cands[-minf] = minx
        (minf, minx) = OneNbin.fit_oi_geom_trunc(
            [0.0, initial],
            f,
            upper = [Inf, 23.0],
            verbose = false
        )
        geom_cands[-minf] = minx
    end
    minx = nb_cands[maximum(keys(nb_cands))]
    N_hat_nb = No + exp(minx[3])
    w_hat_nb = 1.0 / (1.0 + exp(-minx[1]))
    a_hat = exp(minx[2])
    minx = geom_cands[maximum(keys(geom_cands))]
    N_hat_geom = No + exp(minx[2])
    w_hat_geom = 1.0 / (1.0 + exp(-minx[1]))
    out = DataFrame(
        dataset = file,
        No = No,
        N_hat_nb = N_hat_nb,
        w_hat_nb = w_hat_nb,
        a_hat = a_hat,
        N_hat_geom = N_hat_geom,
        w_hat_geom = w_hat_geom
    )
    append!(df, out)
end
CSV.write(OUTPUT_FOLDER * "estimates.csv", df)
