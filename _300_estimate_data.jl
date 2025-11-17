include("_100_utils.jl")
include("_120_one_nbin.jl")
include("_130_benchmarks.jl")

using .OneNbin
using .Utils
using .Benchmarks

using CSV
using DataFrames
using Distributions
using Optim

import Random: seed!

INPUT_FOLDER::String = "./_100_input/datasets/"
OUTPUT_FOLDER::String = "./_900_output/data/datasets/"
SEED::Int = 777
B::Int = 10000 # Number of bootstrap samples for standard error estimation
SIG_LEVEL::Float64 = 0.05 # Significance level for confidence intervals

seed!(SEED)
data_files = [file for file in readdir(INPUT_FOLDER) if endswith(file, ".txt")]
df = DataFrame(
    dataset = String[],
    No = Int64[],
    N_hat_nb = Float64[],
    w_hat_nb = Float64[],
    a_hat = Float64[],
    upper_ci_nb = Float64[],
    lower_ci_nb = Float64[],
    N_hat_geom = Float64[],
    w_hat_geom = Float64[],
    upper_ci_geom = Float64[],
    lower_ci_geom = Float64[]
)
for file in data_files
    println("***$file***")
    file_path = joinpath(INPUT_FOLDER, file)
    new_data = CSV.read(file_path, DataFrame, header=false)
    
    f = Dict(new_data[:, 1] .=> new_data[:, 2])
    No = sum(values(f))
    
    chao = Benchmarks.chao(f)
    initial_N = chao < Inf ? chao : 2 * No

    # Fit MPLE-NB model and store results
    (minf, minx) = OneNbin.fit_oi_nbin_trunc(
        [log(1.0), log(initial_N - No)],
        f;
        verbose = true,
        upper = [10.0, 23.0],
        method = Optim.GradientDescent()
    )

    N_hat_nb = No + exp(minx[2])
    w_hat_nb = OneNbin.w_hat(minx[1], minx[2], f)
    a_hat = exp(minx[1])

    fb = Utils.resample_f(f, N_hat_nb, B)
    Nbs = []
    for b in fb
        res_b = OneNbin.fit_oi_nbin_trunc(
            [log(1.0), log(initial_N - No)],
            b;
            verbose = false,
            upper = [10.0, 23.0],
            method = Optim.GradientDescent()
        )[2]
        push!(Nbs, No + exp(res_b[2]))
    end
    upper_ci_nb = quantile(Nbs, 1.0 - SIG_LEVEL / 2)
    lower_ci_nb = quantile(Nbs, SIG_LEVEL / 2)

    # Fit MPLE-G model and store results
    (minf, minx) = OneNbin.fit_oi_geom_trunc(
        [log(initial_N - No)],
        f,
        verbose = true,
        upper = [23.0],
        method = Optim.LBFGS()
    )

    N_hat_geom = No + exp(minx[1])
    w_hat_geom = OneNbin.w_hat(log(1.0), minx[1], f)

    fb = Utils.resample_f(f, N_hat_geom, B)
    Nbs = []
    for b in fb
        res_b = OneNbin.fit_oi_geom_trunc(
            [log(initial_N - No)],
            b;
            verbose = false,
            upper = [23.0],
            method = Optim.LBFGS()
        )[2]
        push!(Nbs, No + exp(res_b[1]))
    end
    upper_ci_geom = quantile(Nbs, 1.0 - SIG_LEVEL / 2)
    lower_ci_geom = quantile(Nbs, SIG_LEVEL / 2)

    # Store results in DataFrame
    out = DataFrame(
        dataset = file,
        No = No,
        N_hat_nb = N_hat_nb,
        w_hat_nb = w_hat_nb,
        a_hat = a_hat,
        upper_ci_nb = upper_ci_nb,
        lower_ci_nb = lower_ci_nb,
        N_hat_geom = N_hat_geom,
        w_hat_geom = w_hat_geom,
        upper_ci_geom = upper_ci_geom,
        lower_ci_geom = lower_ci_geom
    )
    append!(df, out)
end
CSV.write(OUTPUT_FOLDER * "estimates.csv", df)
