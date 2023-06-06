include("_110_estimator.jl")
include("_120_benchmarks.jl")
include("_100_utils.jl")


import .Estimator
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
n = repeat([5, 100], 20)[1:T]
trials = 500
d = truncated(Beta(a, b), upper = 1.0 / maximum(n))
res = zeros(trials, 15)
ngrid = 75
seed!(7)
p = rand(d, N)
while abs(sum(p) - 1.0) > 0.01
    global p = rand(d, N)
end
output_file = "poisson_estimates.csv"
#Utils.write_row(output_file, ["a_hat", "Nu_hat", "N_hat", "No", "trial", "T", "alpha", "N", "type"])
for i in 1:trials
    S = [sample(1:N, pweights(p * n[t]), n[t], replace = false) for t in 1:T]
    K = Dict{Int, Int}()
    for s in S
        addcounts!(K, s)
    end
    f = countmap(values(K))
    O = Set([i for j in S for i in j])
#    X = Dict{Any, Any}()
#    for i in O
#        X[i] = [i in s for s in S]
#    end
#    println("$(length(X))")
#    (minf, minx, ret) = Estimator.fit_model(X, n, ngrid, [5.0, length(X)]; ftol = 1e-6, upper = [Inf, 10000])
#    mc_estimate = logistic.(rand(truncated(Normal(minx[1], exp(minx[2])), upper = logit(1 / n)), 100000))
#    N_hat = Estimator.u_size(minx[1], (minx[2] + length(X) - 1) * minx[1], maximum(n), 0)
    No = length(O)
    sum_x = collect(values(K))
    theta = [log(5.0), log(length(O))]
    for i in 1:10000
        grad = [Estimator.gradient_a(theta[2], theta[1], n, No, sum_x), Estimator.gradient_Nu(theta[2], theta[1], n, No, sum_x)]
        if sum(isnan.(grad)) > 0
            break
        else
            theta = theta + 0.01 * grad
        end
    end
    N_hat = exp(theta[2]) + No
    draws = []
    push!(draws, [exp(theta[1]), exp(theta[2]), N_hat, length(O), i, T, a, N, "Pseudolikelihood"])
#    row = [exp(theta[1]), N_hat]
#     benchmarks = []
#     push!(benchmarks, Benchmarks.schnabel(S, n))
#     push!(benchmarks, Benchmarks.chao(length(O), f))
#     push!(benchmarks, Benchmarks.zelterman(length(O), f))
#     push!(benchmarks, Benchmarks.conway_maxwell(length(O), f))
#     push!(benchmarks, Benchmarks.huggins(T, K))
#     push!(benchmarks, Benchmarks.turing_geometric(length(O), f, T))
#     push!(benchmarks, Benchmarks.turing(length(O), f, T))
#     for k in 1:5
#         jk = Benchmarks.jackknife(length(O), T, f, k)
#         push!(benchmarks, jk)
#     end
#     for b in benchmarks
#         push!(row, b)
#     end
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
        push!(draws, [-999.0, -999.0, benchmarks[b], length(O), i, T, a, N, b])
    end
    for d in draws
        Utils.write_row(output_file, d)
    end
    #push!(row, length(O))
    #res[i,:] = row
    println(draws[1])
end
