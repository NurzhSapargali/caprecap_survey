module Estimator

import Distributions: Beta, logpdf

using NLopt
using StatsBase

export loglh, fit_model

function log_prior(p::Real, a::Real, b::Real)
    d = Beta(a, b)
    return logpdf(d, p)
end

function datalh(p::Real, x_i::Vector{Bool}, n::Vector{Int})
    comp_prob = 1.0 .- n * p
    return prod((n * p).^x_i .* comp_prob.^(1 .- x_i))
end

function trapezoid(a::Real, b::Real, data::Matrix, grid::LinRange{Float64, Int})
    delta = grid[2] - grid[1]
    prior_vals = [2.0 * exp(log_prior(p, a, b)) for p in grid]
    prior_vals[1] /= 2.0
    prior_vals[lastindex(grid)] /= 2.0
    integrand_vals = transpose(prior_vals) * data
    return transpose(delta / 2.0 * integrand_vals)
end

# function monte_carlo(prob_matrix::Matrix{Real},
#                      complement_prob_matrix::Matrix{Real},
#                      x_i::Vector{Bool})
#     X = repeat(x_i, 1, size(prob_matrix)[2]);
#     G_X = (prob_matrix .^ X) .* (complement_prob_matrix .^ (1 .- X));
#     integrands = prod(G_X, dims=1);
#     return mean(integrands);
# end

function loglh(alpha::Real,
               N_u::Real,
               N_o::Int,
               data::Matrix,
               grid::LinRange{Float64, Int})
    beta = alpha * (N_u + N_o - 1.0)
#     points = rand(Beta(alpha, beta), draws);
#     P = n * transpose(points)
#     comp_P = 1.0 .- P
#     I = keys(X) .|> (g -> monte_carlo(P, comp_P, X[g]))
    cols = size(data)[2]
    I = trapezoid(alpha, beta, data[:, 1:(cols - 1)], grid)
    fails, avg_fail = (length(I[I .< 0]) / length(I), mean(I[I .< 0]))
    I[I .< 0] .= 5e-200
#     truncation = 1.0 - monte_carlo(P, comp_P, zeros(Bool, T));
    unobserved = reshape(data[:, cols], length(data[:, cols]), 1)
    truncation = 1.0 - trapezoid(alpha, beta, unobserved, grid)[1]
    if truncation < 0
        truncation = 1.0
    end
    lh = -N_o * log(truncation) + sum(log.(I))
    println("....alpha = $alpha, N_u = $N_u, error_rate = $fails, avg_error = $avg_fail, lh = $lh")
    return lh;
end

function fit_model(S::Vector,
                   O::Set,
                   n::Vector{Int64},
                   ngrid::Int)
    grid = LinRange(0.0001, 0.9999, ngrid)
    println("Setting up the design matrix....")
    X = Dict{Any, Vector{Bool}}()
    for i in O
        X[i] = [i in s for s in S]
        println("....$(length(O) - length(X)) left")
    end
    N_o = length(X)
    D = reduce(hcat, [[datalh(p, X[i], n) for p in grid] for i in keys(X)])
    D = hcat(D, [datalh(p, zeros(Bool, length(n)), n) for p in grid])
    LL(x, grad) = -loglh(x[1], x[2], N_o, D, grid);
    opt = Opt(:LN_SBPLX, 2);
    lower = [0.01, 0];
    upper = [10000, Inf];
    opt.upper_bounds = upper
    opt.lower_bounds = lower;
    opt.min_objective = LL;
    opt.xtol_abs = 1e-2
    println("Optimizing....")
    (minf, minx, ret) = NLopt.optimize(opt, [5.0, length(O)]);
    return (minf, minx, ret);
end

end
