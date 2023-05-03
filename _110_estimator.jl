module Estimator

import Distributions: Gamma, Beta, logpdf, truncated, logcdf
import StatsFuns: logistic, logit

using NLopt
using StatsBase
using QuadGK

export loglh, log_datalh, log_prior, log_posterior, fit_model, trapezoid, monte_carlo, u_size #, fit_univariate_model

function u_size(a::Real, b::Real, max_n::Int, o_size::Int)
    first = log(a + b) - log(a)
    second = logcdf(Beta(a + 1.0, b), 1.0 / max_n) - logcdf(Beta(a, b), 1.0 / max_n)
    return exp(first + second) - o_size
end

function log_prior(eta::Real, a::Real, b::Real, max_n::Int)
    p = logistic(eta)
#    b = (N - 1.0) * a
    d = truncated(Beta(a, b), upper = 1.0 / max_n)
    return logpdf(d, p) + log(p) + log(1.0 - p)
end

function log_datalh(eta::Real, x_i::Vector{Bool}, n::Vector{Int})
    p = logistic(eta)
    comp_prob = 1.0 .- n * p
    lhs = x_i .* (log.(n) .+ log(p)) .+ (1 .- x_i) .* log.(comp_prob)
    return sum(lhs)
end

function log_posterior(eta::Real, a::Real, b::Real, x_i::Vector{Bool}, n::Vector{Int})
    max_n = maximum(n)
    return log_prior(eta, a, b, max_n) + log_datalh(eta, x_i, n)
end

function trapezoid(a::Real, b::Real, data::Matrix, n::Vector{Int}, grid::LinRange{Float64, Int})
    delta = grid[2] - grid[1]
    prior_vals = log_prior.(grid, a, b, maximum(n))
    I = deepcopy(data)
    for c in 1:size(data)[2]
        I[:,c] += prior_vals
    end
    I[1,:] /= 2.0
    I[size(I)[1],:] /= 2.0
    return sum(2 * exp.(I), dims = 1) * delta / 2.0
end

function monte_carlo(prob_matrix::Matrix,
                     complement_prob_matrix::Matrix,
                     x_i::Vector{Bool})
    X = repeat(x_i, 1, size(prob_matrix)[2]);
    G_X = X .* log.(prob_matrix) .+ (1 .- X) .* log.(complement_prob_matrix);
    integrands = sum(G_X, dims=1);
    return mean(exp.(integrands));
end

function loglh(a::Real,
#               N_u::Real,
               b::Real,
#               X::Dict,
               data::Matrix,
               n::Vector{Int},
#               draws::Int,
               grid::LinRange{Float64, Int},
               verbose::Bool = true)
    try
        N_o = size(data)[2] - 1
#        N_o = length(X)
#        b = a * (N_o + N_u - 1.0)
#        N = N_o + N_u
#        I = [trapezoid(alpha, b, X[i], n, grid) for i in keys(X)]
#        unobserved = trapezoid(alpha, b, zeros(Bool, length(n)), n, grid)
#        beta = alpha * (N_u + N_o - 1.0)
#        z = rand(Gamma(beta), draws)
#        n_max = maximum(n)
#        y = [rand(truncated(Gamma(alpha), upper = exp(logit(1 / n_max)) * i)) for i in z]
#        points = log.(y) - log.(z)
#        d = truncated(Beta(alpha, beta), upper = 1.0 / maximum(n))
#        points = rand(d, draws)
#        P = n * logistic.(transpose(points))
#         if (length(P[P .>= 1]) != 0)
#             P[P .>= 1] .= 0.99
#         end
#        comp_P = 1.0 .- P
#        obs = [monte_carlo(P, comp_P, X[g]) for g in keys(X)]
    #    cols = size(data)[2]
#        upper = logit(1.0 / maximum(n))
#        obs = [quadgk(p -> exp(log_posterior(p, a, b, X[i], n)), -709.99, upper)[1] for i in keys(X)]
#        no_inc = zeros(Bool, length(n))
#        nobs = quadgk(p -> exp(log_posterior(p, a, b, no_inc, n)), -709.99, upper)[1]
        N_u = u_size(a, b, maximum(n), N_o)
#         fails, avg_fail = (length(I[I .< 0]) / length(I), mean(I[I .< 0]))
#         I[I .<= 0] .= 5e-200
#        truncation = 1.0 - monte_carlo(P, comp_P, zeros(Bool, length(n)));
    #    unobserved = reshape(data[:, cols], length(data[:, cols]), 1)
    #    integral_unobserved = 0.0
    #    try
    #        integral_unobserved = quadgk(p -> joint(p, alpha, beta, zeros(Bool, length(n)), n), 0, 1)[1]
    #    catch y
    #        if isa(y, DomainError)
    #            points = rand(Beta(alpha, beta), 100000)
    #            P = n * transpose(points)
    #            comp_P = 1.0 .- P
    #            integral_unobserved = monte_carlo(P, comp_P, zeros(Bool, length(n)))
    #        else
    #            error("Something wrong by truncation")
    #        end
    #    end
        integrals = trapezoid(a, b, data, n, grid)
        obs = integrals[:,1:(lastindex(integrals) - 1)]
        truncation = 1.0 - integrals[1,lastindex(integrals)]
#         bad_truncation = 1.0
#         if truncation <= 0.0
#             bad_truncation = truncation
#             truncation = 5e-200
#         end
         lh = -N_o * log(truncation) + sum(log.(obs))
#        lh = -N_o * log(1.0 - nobs) + sum(log.(obs))
        if verbose
#             if bad_truncation <= 0.0
#                 println("....alpha = $alpha, N_u = $N_u, error_rate = $fails, avg_error = $avg_fail, bad_truncation = $bad_truncation, lh = $lh")
#             else
#                 println("....alpha = $alpha, N_u = $N_u, error_rate = $fails, avg_error = $avg_fail, lh = $lh")
#             end
              println("....alpha = $a, Nu = $N_u, b = $b, lh = $lh")
        end
        return lh
    catch e
        bt = catch_backtrace()
        showerror(stdout, e, bt)
        rethrow(e)
    end
end

function fit_model(X::Dict,
                   n::Vector{Int64},
#                   draws::Int,
                   ngrid::Int,
                   theta0::Vector;
                   lower::Vector = [0.01, 0.01],
                   upper::Vector = [Inf, Inf],
                   ftol::Real = 0.001)
#     grid = LinRange(0.0001, 0.9999, ngrid)
    upper_bound = logit(1.0 / maximum(n) - eps())
    grid = LinRange(-709.99, upper_bound, ngrid)
    D = reduce(hcat, [[log_datalh(eta, X[i], n) for eta in grid] for i in keys(X)])
    D = hcat(D, [log_datalh(eta, zeros(Bool, length(n)), n) for eta in grid])
    LL(x, grad) = -loglh(x[1], x[2], D, n, grid)
    opt = Opt(:LN_SBPLX, 2)
    opt.upper_bounds = upper
    opt.lower_bounds = lower
    opt.min_objective = LL
    opt.ftol_abs = ftol
#    inequality_constraint!(opt, (x,g) -> -u_size(x[1], x[2], maximum(n), length(X)), 0.0)
    println("Optimizing....")
    (minf, minx, ret) = NLopt.optimize(opt, theta0)
    return (minf, minx, ret)
end

# function fit_univariate_model(S::Vector,
#                               O::Set,
#                               n::Vector{Int64},
#                               fixed_val::Real,
#                               ngrid::Int,
#                               fixed_name::String = "alpha")
#     grid = LinRange(0.0001, 0.9999, ngrid)
#     println("Setting up the design matrix....")
#     X = Dict{Any, Vector{Bool}}()
#     for i in O
#         X[i] = [i in s for s in S]
#         println("....$(length(O) - length(X)) left")
#     end
#     N_o = length(X)
#     D = reduce(hcat, [[datalh(p, X[i], n) for p in grid] for i in keys(X)])
#     D = hcat(D, [datalh(p, zeros(Bool, length(n)), n) for p in grid])
#     opt = Opt(:LN_SBPLX, 1)
#     if fixed_name == "alpha"
#         LL1(x, grad) = -loglh(fixed_val, x[1], N_o, D, grid);
#         lower = [0]
#         upper = [Inf]
#         initial_guess = [N_o]
#         opt.min_objective = LL1
#     elseif fixed_name == "Nu"
#         LL2(x, grad) = -loglh(x[1], fixed_val, N_o, D, grid);
#         upper = [10000]
#         lower = [0.01]
#         initial_guess = [5.0]
#         opt.min_objective = LL2
#     end
#     opt.upper_bounds = upper
#     opt.lower_bounds = lower
#     opt.xtol_abs = 1e-2
#     println("Optimizing....")
#     (minf, minx, ret) = NLopt.optimize(opt, initial_guess)
#     return (minf, minx, ret)
# end

end
