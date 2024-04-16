include("./_100_utils.jl")
include("./_120_gamma_estimator.jl")
include("./_130_benchmarks.jl")

import .GammaEstimator
import .Benchmarks
import .Utils

using DataFrames
using CSV

const DATA_FOLDER = "./_900_output/data/user_nets/"
const INDICES = [250, 251, 252, 253, 254, 980, 981, 982, 985, 986]

function read_user_net(file)
    df = DataFrame(CSV.File(file))
    #df = df[df[:, :i] .!= df[:, :j], :]
    return Set(vcat(df[:, :i], df[:, :j]))
end

nov_2020 = []
nov_2022 = []
for ind in INDICES
    file = DATA_FOLDER * "user_net_$ind.csv"
    capture = read_user_net(file)
    if ind > 900
        push!(nov_2022, capture)
    else
        push!(nov_2020, capture)
    end
end

estimates = DataFrame(
    "a_hat" => Float64[],
    "pseudo" => Float64[],
    "schnabel" => Float64[],
    "chao" => Float64[],
    "zelterman" => Float64[],
    "cmp" => Float64[],
    "turing_geometric" => Float64[],
    "turing" => Float64[],
    "morgan_ridout" => Float64[],
    "chao_lee_jeng_0" => Float64[],
    "chao_lee_jeng_1" => Float64[],
    "chao_lee_jeng_2" => Float64[],
    "jackknife_1" => Float64[],
    "jackknife_2" => Float64[],
    "jackknife_3" => Float64[],
    "jackknife_4" => Float64[],
    "jackknife_5" => Float64[],
    "observed" => Int[]
)

for S in [nov_2020, nov_2022]
    T = length(S)
    n = [length(s) for s in S]
    f = Utils.freq_of_freq(Utils.cap_freq(S))
    No = sum(values(f))
    (minf, minx, ret) = GammaEstimator.fit_Gamma(
        [log(1.0), log(1.0)],
        n,
        No,
        f,
        ftol = 1e-15,
        lower = [log(0.0001), -Inf],
        upper = [10.0, 30]
    )
    row = [
        exp(minx[1]),
        exp(minx[2]) + sum(values(f)),
        Benchmarks.schnabel(S, n),
        Benchmarks.chao(No, f),
        Benchmarks.zelterman(No, f),
        Benchmarks.conway_maxwell(No, f),
        Benchmarks.turing_geometric(No, f, T),
        Benchmarks.turing(No, f, T),
        Benchmarks.morgan_ridout(f, T, "./estimateN.R")
    ]
    for b in 0:2
        push!(row, Benchmarks.chao_lee_jeng(No, f, T, n, b))
    end
    for k in 1:5
        jk = Benchmarks.jackknife(No, T, f, k)
        push!(row, jk)
    end
    push!(row, sum(values(f)))
    push!(estimates, row)
end

estimates[!, :date] = ["Nov 2020", "Nov 2022"]
CSV.write("./_900_output/data/user_nets/estimates.csv", estimates)