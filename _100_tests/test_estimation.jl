include("../_120_gamma_estimator.jl")
include("../_110_beta_estimator.jl")
include("../_130_benchmarks.jl")
include("../_100_utils.jl")

import .GammaEstimator
import .BetaEstimator
import .Benchmarks
import .Utils

using Distributions
using StatsBase
using Plots
using StatsFuns

import Random: seed!

N = 10000
a = 1.0
b = a * (N - 1.0)
T = 20
n = repeat([37, 2, 100, 17, 2, 75, 17, 44, 112, 3], 5)[1:T]
trials = 1
d = truncated(Beta(a, b), upper = 1.0 / maximum(n))
res = zeros(trials, 18)
ngrid = 75
seed!(7)
p = rand(d, N)
S = [Utils.ar_pareto_sample(p, n[t]) for t in 1:T]
K = Dict{Int, Int}()
for s in S
    addcounts!(K, s)
end
f = countmap(values(K))
O = Set([i for j in S for i in j])
X = Dict{Any, Any}()
for i in O
    X[i] = [i in s for s in S]
end
println("$(length(X))")
x_sums = Dict(i => sum(X[i]) for i in keys(X))
(minf1, minx1, ret1) = GammaEstimator.fit_Gamma(
    [5.0, Benchmarks.turing(length(O), f, T) - length(X)],
    n,
    length(X),
    x_sums
)
N_hat1 = length(X) + minx1[2]
row = [N_hat1]
(minf2, minx2, ret2) = BetaEstimator.fit_Beta(X, n, ngrid, [5.0, Benchmarks.turing(length(O), f, T) - length(X)]; ftol = 1e-4, upper = [200, Inf])
N_hat2 = minx2[2] + length(X)
push!(row, N_hat2)
benchmarks = Dict{}()
benchmarks["Schnabel"] = Benchmarks.schnabel(S, n)
benchmarks["Chao"] = Benchmarks.chao(length(O), f)
benchmarks["Zelterman"] = Benchmarks.zelterman(length(O), f)
benchmarks["Conway-Maxwell-Poisson"] = Benchmarks.conway_maxwell(length(O), f)
benchmarks["Turing Geometric"] = Benchmarks.turing_geometric(length(O), f, T)
benchmarks["Turing"] = Benchmarks.turing(length(O), f, T)
MR_hats = Benchmarks.morgan_ridout(f, T, "../estimateN.R")
benchmarks["Morgan Ridout"] = MR_hats[argmin(abs.(MR_hats .- N))]
for b in 0:2
    benchmarks["Chao Lee Jeng $b"] = Benchmarks.chao_lee_jeng(length(O), f, T, n, b)
end
for k in 1:5
    jk = Benchmarks.jackknife(length(O), T, f, k)
    benchmarks["Jackknife k = $(k)"] = jk
end
for b in keys(benchmarks)
    push!(row, benchmarks[b])
end
push!(row, length(X))
res[1,:] = row
println(res[1,:])
