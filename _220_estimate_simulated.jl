include("_100_utils.jl")
include("_120_one_nbin.jl")
include("_130_benchmarks.jl")

import .OneNbin
import .Utils
import .Benchmarks

using StatsBase

import Random: seed!

ALPHAS::Vector{Float64} = [2.0]#[0.5, 2.0, 10.0]
DATA_FOLDER::String = "./_200_input/simulated/"
breaks_T::Vector{Int64} = collect(5:5:50)
OUTPUT_FOLDER::String = "./_900_output/data/simulated/"
pops::Vector{Int64} = [1000, 5000, 10000]
SEED::Int = 777

seed!(SEED)
for i in eachindex(ALPHAS)
    alpha = ALPHAS[i]
    output_file = OUTPUT_FOLDER * "estimates_$(alpha).csv"
    #Utils.write_row(output_file, ["w_hat", "a_hat", "Nu_hat", "N_hat", "No", "trial", "T", "alpha", "N", "type"])
    for N in pops
        data_folder = DATA_FOLDER * "alpha_$(alpha)/"
        #data_files = [file for file in readdir(data_folder) if occursin("sample", file) && occursin("$(N).csv", file)]
        data_files = ["sample_256_$(N).csv", "sample_331_$(N).csv"]
        for file in data_files
            trial_no = parse(Int, split(split(file, "_")[2], ".")[1])
            samples = Utils.read_captures(data_folder * file)
            draws = Array{Any}(undef, 160)
            Threads.@threads for i in eachindex(breaks_T)
                j = 1
                t = breaks_T[i]
                S = samples[1:t]
                K = Utils.cap_freq(S)
                f = Utils.freq_of_freq(K)
                println("***TRIAL NO $file, $t***")
                n = [length(s) for s in S]
                No = sum(values(f))
                benchmarks = Dict{}()
                benchmarks["Turing"] = Benchmarks.turing(No, f, t)
                (minf, minx) = OneNbin.fit_oi_nbin_trunc(
                    [0.0, log(1.0), log(benchmarks["Turing"] - No)],
                    f,
                    upper = [Inf, 20.0, 23.0]
                )
                N_hat = No + exp(minx[3])
                w = 1.0 / (1.0 + exp(-minx[1]))
                println([minx[1], minx[2], minx[3], N_hat, No, trial_no, t, alpha, N])
                draws[(i - 1) * 16 + j] = [w, exp(minx[2]), exp(minx[3]), N_hat, No, trial_no, t, alpha, N, "MPLE-NB"]
                j += 1
                (minf, minx) = OneNbin.fit_oi_geom_trunc(
                    [0.0, log(benchmarks["Turing"] - No)],
                    f,
                    upper = [Inf, 23.0]
                )
                N_hat = No + exp(minx[2])
                w = 1.0 / (1.0 + exp(-minx[1]))
                draws[(i - 1) * 16 + j] = [w, 1.0, exp(minx[2]), N_hat, No, trial_no, t, alpha, N, "MPLE-G"]
                j += 1
                benchmarks["Schnabel"] = Benchmarks.schnabel(S, n)
                benchmarks["Chao"] = Benchmarks.chao(No, f)
                benchmarks["Zelterman"] = Benchmarks.zelterman(No, f)
                benchmarks["Conway-Maxwell-Poisson"] = Benchmarks.conway_maxwell(No, f)
                benchmarks["Turing Geometric"] = Benchmarks.turing_geometric(No, f, t)
                for b in 0:2
                    benchmarks["Chao Lee Jeng $b"] = Benchmarks.chao_lee_jeng(No, f, t, n, b)
                end
                for k in 1:5
                    jk = Benchmarks.jackknife(No, t, f, k)
                    benchmarks["Jackknife k = $(k)"] = jk
                end
                for b in keys(benchmarks)
                    draws[(i - 1) * 16 + j] = [-999.0, -999.0, benchmarks[b] - No, benchmarks[b], No, trial_no, t, alpha, N, b]
                    j += 1
                end
            end
            for t in breaks_T
                S = samples[1:t]
                K = Utils.cap_freq(S)
                f = Utils.freq_of_freq(K)
                n = [length(s) for s in S]
                No = sum(values(f))
                mr_hat = Benchmarks.morgan_ridout(f, t, "./estimateN.R")
                push!(draws, [-999.0, -999.0, mr_hat - No, mr_hat, No, trial_no, t, alpha, N, "Morgan-Ridout"])
            end
            for d in draws
                Utils.write_row(output_file, d)
            end
        end
    end
end
