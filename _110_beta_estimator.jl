module BetaEstimator

import Distributions: Beta, quantile

using Optim
using StatsBase

export approx_loglike, fit_Beta

function cond_likelihood(mu::Vector, x::Vector)
    """
        cond_likelihood(mu, x)

    Compute conditional likelihood of capture history x given capture probability mu
    """
    return prod(mu .^ x .* (1.0 .- mu) .^ (1 .- x))
end

function mc_integration(
    p::Vector,
    x::Vector,
    n::Vector{Int},
)
    """
        mc_integration(p, x, n)

    Compute Monte Carlo integration of likelihood over capture history x
    """
    L(roll) = cond_likelihood(roll * n, x)
    return mean(L.(p))
end

function approx_loglike(
    a::Real,
    Nu::Real,
    X::Matrix,
    u::Vector;
    verbose::Bool = true
)
    """
        loglh(a, N_u, X, u, [verbose,])

    Compute marginal likelihood over entire data at parameter a, number of
    unobserved individuals Nu, capture history of observed individuals X
    and quantiles u for Monte Carlo integration
    """
    No = size(X)[1]
    N = Nu + No
    n = vec(sum(X, dims = 1))

    d = Beta(a, a * (N - 1.0))
    p = quantile.(d, u)
    p = p[p .<= 1.0 / maximum(n)]

    marginal(x) = mc_integration(p, x, n)
    lh = sum(log.(mapslices(marginal, X, dims = 2)))
    no_cap = zeros(Int, length(n))
    lh += -No * log(1.0 - marginal(no_cap))

    if verbose
        println("a = $a, N = $N, lh = $lh")
    end
    return lh
end

function fit_Beta(
    theta0::Vector,
    X::Matrix;
    lower::Vector = [-Inf, -Inf],
    upper::Vector = [Inf, Inf],
    draws::Int = 10000,
    ftol = 1e-5,
    verbose::Bool = true,
)
    """
        fit_Beta(X, n, ngrid, theta0, [lower, upper, ftol,])

        Maximize likelihood with respect to model parameters
    """
    u = rand(draws)

    L(x) = -approx_loglike(exp(x[1]), exp(x[2]), X, u; verbose = verbose)

    res = optimize(
        L,
        lower,
        upper,
        theta0,
        Fminbox(NelderMead()),
        Optim.Options(iterations = 50, f_tol = ftol)
    )
    (minf, minx) = (Optim.minimum(res), Optim.minimizer(res))
    return (minf, minx)
end

end
