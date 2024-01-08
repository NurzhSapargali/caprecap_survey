module BetaEstimator

import Distributions: Beta
import SpecialFunctions: logbeta

using LinearAlgebra
using NLopt
using StatsBase

export approx_loglh, fit_Beta

function logitBeta(y, a::Real, b::Real, c::Real)
    p = 1.0 / (exp(-y) + 1.0)
    logfp = (a - 1.0) * log(p) + (b - 1.0) * log(c - p) - (a + b - 1.0) * log(c) - logbeta(a, b)
    logdy = log(p) + log(1.0 - p)
    return exp(logfp + logdy)
end

function joint_density(y, x::Vector, a::Real, b::Real, n::Vector)
    p = 1.0 / (exp(-y) + 1.0)
    c = 1.0 / maximum(n)
    loglike = sum(x .* log.(n * p) + (1.0 .- x) .* log.(1.0 .- n * p))
    dsty = log(logitBeta.(y, a, b, c)) + loglike
    return exp(dsty)
end

function approx_loglh(
    a::Real,
    Nu::Real,
#    X::Matrix,
    X::Dict,
    n::Vector{Int},
    grid_size::Int;
    verbose::Bool = true
)
    """
        loglh(a, N_u, X, n, ngrid, [verbose,])

    Compute marginal likelihood over entire data at parameter a, number of
    unobserved individuals N_u, capture history of observed individuals X,
    sample sizes n and length of grid of quantiles for numerical integration
    """
    #No = size(X)[1] - 1
    No = length(X)
    N = Nu + No
    n_max_inv = 1.0 / maximum(n)
    b = (N * n_max_inv - 1.0) * a
    dist = Beta(a, b)
    qr = LinRange(0.0001, 0.9999, grid_size)
    #grid_points = quantile(d, qr) * n_max_inv
    logit_grid = log.(
        (quantile(dist, qr) * n_max_inv) ./ (1.0 .- (quantile(dist, qr) * n_max_inv))
    )
    # Perform check for gaps
    dsty = logitBeta.(logit_grid, a, b, n_max_inv)
    logit_grid = logit_grid[.!isnan.(dsty)]
    dsty = dsty[.!isnan.(dsty)]
    remain = lastindex(logit_grid)
    del = logit_grid[2:remain] - logit_grid[1:remain - 1]
    prior_int = sum((dsty[2:remain] + dsty[1:remain - 1]) .* del / 2.0)
    logsum = 0
    for i in keys(X)
        x = get(X, i, repeat([0], length(n)))
        integrands = [joint_density(j, x, a, b, n) for j in logit_grid]
        integral = sum((integrands[1:remain - 1] + integrands[2:remain]) .* del / 2.0)
        #println(integral)
        logsum += log(integral)
    end
    x = repeat([0], length(n))
    integrands = [joint_density(j, x, a, b, n) for j in logit_grid]
    truncation = sum((integrands[1:remain - 1] + integrands[2:remain]) .* del / 2.0)
    #D = transpose(repeat(points, 1, length(n)))
    #P = diagm(n) * D
    #P[P .== 1.0] .= 0.999999999999999
    #P[P .== 0.0] .= 1e-323
    #integrals = sum(exp.(X * log.(P) + (1 .- X) * log.(1.0 .- P)), dims = 2) / size(D)[2]
    lh = -No * log(1.0 - truncation) + logsum
    if verbose
        println("a = $a, b = $b, N = $N, prior integral = $prior_int, lh = $lh")
    end
    return lh
end

function fit_Beta(
    theta0::Vector,
    #X::Matrix,
    X::Dict,
    n::Vector{Int},
    grid_size::Int;
    lower::Vector = [0.01, 0.1],
    upper::Vector = [1000, 100000],
    xtol = 1e-3
)
    """
        fit_model(X, n, ngrid, theta0, [lower, upper, ftol,])

        Maximize likelihood with respect to model parameters
    """

    LL(x, grad) = -approx_loglh(x[1], x[2], X, n, grid_size)
    opt = Opt(:LN_SBPLX, 2)
    opt.min_objective = LL
    opt.xtol_rel = xtol
    opt.lower_bounds = lower
    opt.upper_bounds = upper
    opt.maxeval = 10000
    #opt.vector_storage = 10
    (minf,minx,ret) = optimize(opt, theta0)
    return (minf, minx, ret)
end

end
