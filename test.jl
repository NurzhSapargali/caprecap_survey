include("_110_estimator.jl")
# include("_100_utils.jl")
include("_120_benchmarks.jl")
# 
# using .Estimator
# using .Utils
# using .Benchmarks
# 
# using StatsBase
# using DelimitedFiles
# 
# import Random: seed!
# import Distributions: Gamma, Beta, logpdf, logcdf, truncated
# import StatsFuns: logistic, logit
# import Plots: plot, histogram, plot!, vline!
# 
# ALPHAS::Vector{Float64} = [0.5, 1.0, 5.0, 10.0]
# DATA_FOLDER::String = "./_200_input/diffp/"
# breaks_T::Vector{Int64} = [5, 10, 15, 20]
# OUTPUT_FOLDER::String = "./_900_output/data/diffp/"
# 
# 
# function get_truth(filename)
#     data = readdlm(filename, ' ', String, '\n')[2, 1]
#     return parse(Float64, split(data, ",")[1])
# end
# 
# function logN(a, b, upper)
#     main_term = (log(a + b) - log(a))
#     temper = logcdf(Beta(a, b), upper) - logcdf(Beta(a + 1, b), upper)
#     println(exp(main_term))
#     println(exp(temper))
#     return main_term + temper
# end
# 
# seed!(100)
# alpha = 0.5
# data_folder = DATA_FOLDER * "alpha_$(alpha)/"
# data_files = [file for file in readdir(data_folder) if occursin("perfect_sample", file)]
# N = get_truth(data_folder * "metadata_$(alpha).csv")
# file = "sample_10.csv"
# samples = read_captures(data_folder * file)
# t = 20
# S = samples[1:t]
# K = Dict{Int, Int}()
# for s in S
#     addcounts!(K, s)
# end
# f = countmap(values(K))
# println("***TRIAL NO $file, $t***")
# O = Set([i for j in S for i in j])
# n = [length(s) for s in S]
# X = Dict{Any, Vector{Bool}}()
# for i in O
#     X[i] = [i in s for s in S]
#     println("....$(length(O) - length(X)) left")
# end
# upper_bound = logit(1 / maximum(n)) - eps()
# grid = LinRange(-700, upper_bound, 15000)
# beta = alpha * (N - 1.0)
# z = rand(Gamma(beta), 3500)
# n_max = maximum(n)
# y = [rand(truncated(Gamma(alpha), upper = exp(logit(1 / n_max)) * i)) for i in z]
# points = log.(y) - log.(z)
# histogram(points, normalize = :pdf, label = "")
# plot!(grid, exp.(log_prior.(grid, 0.5, N, maximum(n))), xlim = [-20, Inf], label = "")
# vline!([logit(1 / maximum(n)), -log(N - 1)], label = "")
# P = n * logistic.(transpose(points))
# comp_P = 1.0 .- P
# mc = [monte_carlo(P, comp_P, X[g]) for g in keys(X)]
# D = reduce(hcat, [[log_datalh(eta, X[i], n) for eta in grid] for i in keys(X)])
# quad = transpose(trapezoid(0.5, N, D, n, grid))
# histogram(quad - mc, title = "Trapezoid (15k grid points) - MC (3.5k draws)", label = "")
# (minf, minx, ret) = fit_model(X, n, [5.0, length(X)]; ftol = 0.000000001)
# println(exp(logN(minx[1], minx[2], 1 / maximum(n))))

import .Estimator
import .Benchmarks

using Distributions
using StatsBase
using Plots
using StatsFuns

import Random: seed!

N = 1000
a = 0.5
b = a * (N - 1.0)
T = 5
n = repeat([100, 150], 10)[1:T]
trials = 1
d = truncated(Beta(a, b), upper = 1.0 / maximum(n))
res = zeros(trials, 13)
ngrid = 75
seed!(7)
for i in 1:trials
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
    (minf, minx, ret) = Estimator.fit_model(X, n, ngrid, [5.0, length(X)]; ftol = 1e-4, upper = [Inf, Inf])
#    mc_estimate = logistic.(rand(truncated(Normal(minx[1], exp(minx[2])), upper = logit(1 / n)), 100000))
    N_hat = Estimator.u_size(minx[1], (minx[2] + length(X) - 1) * minx[1], maximum(n), 0)
    row = [N_hat]
    benchmarks = Dict{}()
    benchmarks["Schnabel"] = Benchmarks.schnabel(S, n)
    benchmarks["Chao"] = Benchmarks.chao(length(O), f)
    benchmarks["Zelterman"] = Benchmarks.zelterman(length(O), f)
    #benchmarks["Conway-Maxwell-Poisson"] = Benchmarks.conway_maxwell(length(O), f)
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
    res[i,:] = row
    println(res[i,:])
end
