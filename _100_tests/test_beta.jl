include("../_110_beta_estimator.jl")
include("../_130_benchmarks.jl")

import .BetaEstimator
import .Benchmarks

using Distributions
using StatsBase
using Plots
using StatsFuns

import Random: seed!

N = 10000
a = 1.0
b = a * (N - 1.0)
T = 15
n = repeat([37, 2, 202, 17, 2, 75, 17, 44, 112, 3], 2)[1:T]
trials = 1
d = truncated(Beta(a, b), upper = 1.0 / maximum(n))
res = zeros(trials, 14)
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
(minf, minx, ret) = BetaEstimator.fit_Beta(X, n, ngrid, [5.0, length(X)]; ftol = 1e-4, upper = [Inf, Inf])
N_hat = BetaEstimator.u_size(minx[1], (minx[2] + length(X) - 1) * minx[1], maximum(n), 0) + length(X)
row = [N_hat]
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
print(var(values(f)))
