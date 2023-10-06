include("../_120_gamma_estimator.jl")
include("../_110_beta_estimator.jl")
include("../_130_benchmarks.jl")

import .GammaEstimator
import .BetaEstimator
import .Benchmarks

using Distributions
using StatsBase
using Plots
using StatsFuns

import Random: seed!

N = 1000
a = 5.0
b = a * (N - 1.0)
T = 20
n = repeat(Int.(ceil.([37, 2, 100, 17, 2, 75, 17, 44, 112, 3] ./ 1.0)), 2)[1:T]
trials = 1
d = truncated(Beta(a, b), upper = 1.0 / maximum(n))
res = zeros(trials, 15)
ngrid = 75
seed!(7)
p = rand(d, N)
S = [sample(1:N, pweights(p * n[t]), n[t], replace = false) for t in 1:T]
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
sum_x = [sum(i) for i in values(X)]
(minf1, minx1, ret1) = GammaEstimator.fit_Gamma(
    [log(100.0), log(5.0)],
    n,
    length(X),
    sum_x
)
N_hat1 = length(X) + exp(minx1[1])
row = [N_hat1]
(minf2, minx2, ret2) = BetaEstimator.fit_Beta(X, n, ngrid, [5.0, length(X)]; ftol = 1e-4, upper = [Inf, Inf])
N_hat2 = minx2[2] + length(X)
push!(row, N_hat2)
benchmarks = Dict{}()
benchmarks["Schnabel"] = Benchmarks.schnabel(S, n)
benchmarks["Chao"] = Benchmarks.chao(length(O), f)
benchmarks["Zelterman"] = Benchmarks.zelterman(length(O), f)
benchmarks["Conway-Maxwell-Poisson"] = Benchmarks.conway_maxwell(length(O), f)
benchmarks["Huggins"] = Benchmarks.huggins(T, K)
benchmarks["Turing Geometric"] = Benchmarks.turing_geometric(length(O), f, T)
benchmarks["Turing"] = Benchmarks.turing(length(O), f, T)
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
