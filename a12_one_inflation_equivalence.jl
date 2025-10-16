include("_100_utils.jl")
include("_120_one_nbin.jl")
include("_130_benchmarks.jl")

using .OneNbin
using .Utils

using CSV
using DataFrames
using Distributions
using Optim
using Plots

import Random: seed!
import .Benchmarks: chao

INPUT_FOLDER::String = "./_100_input/datasets/"
OUTPUT_FOLDER::String = "./_900_output/data/appendix/one_inflation_equiv/"
UNSTABLE::Vector{String} = ["meth_usage.txt", "pleiades.txt", "illegal_immigrants.txt", "domestic_violence.txt"]
FIGURE_FOLDER::String = "./_900_output/figures/appendix/"


data_files = [file for file in readdir(INPUT_FOLDER) if endswith(file, ".txt")]
df = DataFrame(
    dataset = String[],
    No = Float64[],
    a_hat_oizt = Float64[],
    a_hat_ztoi = Float64[],
    log_pseudo_oizt = Float64[],
    log_pseudo_ztoi = Float64[],
    N_hat_oizt = Float64[],
    N_hat_ztoi = Float64[],
    w_hat_oizt = Float64[],
    w_hat_ztoi = Float64[]
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

    N_hat_oizt = No + exp(minx[2])
    w_hat_oizt = OneNbin.w_hat(minx[1], minx[2], f)
    a_hat_oizt = exp(minx[1])
    log_pseudo_oizt = -minf

    # Fit MPLE-ZTOI model and store results
    (minf, minx) = OneNbin.fit_trunc_nbin_oi(
        [0.5, log(1.0), log(initial_N - No)],
        f,
        verbose = true,
        method = Optim.LBFGS()
    )

    N_hat_ztoi = No + exp(minx[3])
    w_hat_ztoi = minx[1]
    a_hat_ztoi = exp(minx[2])
    log_pseudo_ztoi = -minf

    # Store results in DataFrame
    out = DataFrame(
        dataset = file,
        No = No,
        a_hat_oizt = a_hat_oizt,
        a_hat_ztoi = a_hat_ztoi,
        log_pseudo_oizt = log_pseudo_oizt,
        log_pseudo_ztoi = log_pseudo_ztoi,
        N_hat_oizt = N_hat_oizt,
        N_hat_ztoi = N_hat_ztoi,
        w_hat_oizt = w_hat_oizt,
        w_hat_ztoi = w_hat_ztoi
    )
    append!(df, out)
end
CSV.write(OUTPUT_FOLDER * "one_inf_equiv.csv", df)

cut = df[.!in.(df.dataset, Ref(UNSTABLE)), :]
cut = sort!(cut, [:w_hat_oizt])
plt = plot(
    xlabel = "MPLE of w (OIZT)",
    ylabel = "MPLE of w (ZTOI)",
    size = (1280 / 1.5, 720 / 1.6),
    leftmargin = 5Plots.mm,
    bottommargin = 10Plots.mm,
    xtickfontsize = 12,
    ytickfontsize = 12,
    xguidefontsize = 12,
    yguidefontsize = 12,
    legendfontsize = 12
)
plot!(
    plt,
    cut.w_hat_oizt,
    cut.w_hat_ztoi,
    label = "",
    seriestype = :scatter,
    markersize = 6,
    color = "#003f5c"
)
plot!(
    plt,
    cut.w_hat_oizt,
    cut.w_hat_ztoi,
    label = "",
    linewidth = 2,
    color = "#003f5c",
    linestyle = :dot
)
savefig(
    plt,
    FIGURE_FOLDER * "one_inflation_equiv.pdf"
)