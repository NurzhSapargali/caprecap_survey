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
a = 0.5
T = 200
n = repeat([37, 44, 100, 17, 2, 75, 17, 2, 112, 3], 100)[1:T]
b = a * (N / maximum(n) - 1.0)
trials = 1
d = Beta(a, b)
res = zeros(trials, 17)
ngrid = 75
seed!(7)
p = rand(d, N) / maximum(n)
S = [Utils.ar_pareto_sample(p, n[t]) for t in 1:T]
n = [length(s) for s in S]
T = length(n)
K = Dict{Int, Int}()
for s in S
    addcounts!(K, s)
end
f = Dict(5 => 6, 4 => 21, 6 => 3, 2 => 243, 3 => 58, 1 => 1253) # At T = 50 causes error
f = countmap(values(K))
O = Set([i for j in S for i in j])
converged = false
theta = [log(1.0), log(1.0)]
(minf1, minx1, ret1) = GammaEstimator.fit_Gamma(
    theta,
    f,
    ftol = 1e-7,
    upper = [Inf, Inf],
    lower = [-Inf, -Inf],
)
N_hat1 = length(K) + exp(minx1[2])
# N_hat2 = minx2[2] + length(X)
# push!(row, N_hat2)
benchmarks = Dict{}()
benchmarks["Schnabel"] = Benchmarks.schnabel(S, n)
benchmarks["Chao"] = Benchmarks.chao(sum(values(f)), f)
benchmarks["Zelterman"] = Benchmarks.zelterman(sum(values(f)), f)
benchmarks["Conway-Maxwell-Poisson"] = Benchmarks.conway_maxwell(sum(values(f)), f)
benchmarks["Turing Geometric"] = Benchmarks.turing_geometric(sum(values(f)), f, T)
benchmarks["Turing"] = Benchmarks.turing(sum(values(f)), f, T)
benchmarks["Morgan Ridout"] = Benchmarks.morgan_ridout(f, T, "./estimateN.R")
for b in 0:2
    benchmarks["Chao Lee Jeng $b"] = Benchmarks.chao_lee_jeng(length(O), f, T, n, b)
end
for k in 1:5
    jk = Benchmarks.jackknife(sum(values(f)), maximum(keys(f)), f, k)
    benchmarks["Jackknife k = $(k)"] = jk
end
println(benchmarks)
row = [N_hat1]
for b in keys(benchmarks)
    push!(row, benchmarks[b])
end
push!(row, length(K))
res[1,:] = row
println(res[1,:])
#check = BetaEstimator.fit_Beta(
#    [5.0, benchmarks["Turing"] - length(X)],
#    X,
#    n,
#   100
#)
