module BetaEstimator

import Distributions: Beta, logpdf, truncated, logcdf
import StatsFuns: logistic, logit

using NLopt
using StatsBase

export loglh, log_datalh, log_prior, log_posterior, fit_model, marginal, u_size

function u_size(a::Real, b::Real, max_n::Int, o_size::Int)
    """
        u_size(a, b, max_n, o_size)

    Calculate estimate of size of unobserved individuals based on prior parameters
    a, b and truncation point 1 / max_n by estimation population size and
    subtracting size of observed individuals
    """
    first = log(a + b) - log(a)
    second = logcdf(Beta(a + 1.0, b), 1.0 / max_n) - logcdf(Beta(a, b), 1.0 / max_n)
    return exp(first + second) - o_size
end

function log_prior(eta::Real, a::Real, b::Real, max_n::Int)
    """
        log_prior(eta, a, b, max_n)

    Compute logarithm of logit-beta density with parameters a and b at point
    eta truncated at 1.0 / max_n from above
    """
    p = logistic(eta)
    d = truncated(Beta(a, b), upper = 1.0 / max_n)
    return logpdf(d, p) + log(p) + log(1.0 - p)
end

function small_grid(bqs::Vector{Float64}, a::Real, b::Real, max_n::Int, tol::Real)
    """
        small_grid(bqs, a, b, max_n, tol)

    Perform trapezoid numerical integration over the range of quantiles bqs of
    logit beta distribution with parameters a and b truncated at max_n from above.
    Additionally, check if the accuracy of the numerical integration is within
    tol level
    """
    fun = exp.(log_prior.(logit.(bqs), a, b, max_n))
    deltas = logit.(bqs[2:lastindex(bqs)]) .- logit.(bqs[1:lastindex(bqs) - 1])
    # Compute integral
    I = sum((fun[1:lastindex(bqs) - 1] .+ fun[2:lastindex(bqs)]) / 2.0 .* deltas)
    # Check if accurate
    if abs(1.0 - I) > tol
        return (true, I)
    else
        return (false, I)
    end
end

function log_datalh(eta::Real, x_i::Vector{Bool}, n::Vector{Int})
    """
        log_datalh(eta, x_i, n)

    Compute logarithm of probability mass function of independent non-identically distributed Bernoulli variables with at points of x_i parameters corresponding
    to each element of element-wise product of eta and n
    """
    p = logistic(eta)
    comp_prob = 1.0 .- n * p
    lhs = x_i .* (log.(n) .+ log(p)) .+ (1 .- x_i) .* log.(comp_prob)
    return sum(lhs)
end

function log_posterior(
    eta::Real,
    a::Real,
    b::Real,
    x_i::Vector{Bool},
    n::Vector{Int}
)
    """
        log_posterior(eta, a, b, x_i, n)

    Calculate logarithm of joint density of beta distributed parameters and
    independent non-identically distributed Bernoulli data.
    """
    max_n = maximum(n)
    return log_prior(eta, a, b, max_n) + log_datalh(eta, x_i, n)
end

function marginal(
    a::Real,
    b::Real,
    x_i::Vector{Bool},
    n::Vector{Int},
    bqs::Vector{Float64}
)
    """
        marginal(a, b, x_i, n, bqs)

    Compute marginal density at x_i using trapezoid integration from logarithm
    of posterior desnsity from quantiles of prior bqs and parameters a, b, n
    """
    deltas = logit.(bqs[2:lastindex(bqs)]) - logit.(bqs[1:lastindex(bqs) - 1])
    f = exp.([log_posterior(logit(q), a, b, x_i, n) for q in bqs])
    return sum((f[1:lastindex(f) - 1] .+ f[2:lastindex(f)]) / 2.0 .* deltas)
end

function loglh(
    a::Real,
    N_u::Real,
    X::Dict,
    n::Vector{Int},
    ngrid::Int;
    verbose::Bool = true
)
    """
        loglh(a, N_u, X, n, ngrid, [verbose,])

    Compute marginal likelihood over entire data at parameter a, number of
    unobserved individuals N_u, capture history of observed individuals X,
    sample sizes n and length of grid of quantiles for numerical integration
    """
    try
        N_o = length(X)
        N = N_u + N_o
        b = (N - 1.0) * a
        # Prior quantiles for numerical integration
        start_grid = 0.01
        end_grid = 0.99
        qgrid = LinRange(start_grid, end_grid, ngrid)
        # Truncation upper limit
        prior_dist = truncated(Beta(a, b), upper = 1.0 / maximum(n))
        qs = [quantile(prior_dist, q) for q in qgrid]
        edge_case = false
        # Check if numerical integration is within 0.01 for prior density
        (bad_accuracy, prior_integral) = small_grid(qs, a, b, maximum(n), 0.01)
        while bad_accuracy
            # Keep increasing quantile range over and over
            start_grid = start_grid * 0.1
            end_grid = (1.0 - end_grid) * 0.9 + end_grid
            qgrid = LinRange(start_grid, end_grid, ngrid)
            qs = [quantile(prior_dist, q) for q in qgrid]
            # If the range is too large that it includes 0.0, retreat to default
            if length(qs[qs .== 0.0]) != 0
                obs = 0.0
                truncation = Inf
                N = N_o
                edge_case = true
                break
            end
            (bad_accuracy, prior_integral) = small_grid(qs, a, b, maximum(n), 0.01)
        end
        if !edge_case
            obs = [marginal(a, b, X[i], n, qs) for i in keys(X)]
            truncation = 1.0 - marginal(a, b, zeros(Bool, length(n)), n, qs)
            N = u_size(a, b, maximum(n), 0)
        end
            lh = -N_o * log(truncation) + sum(log.(obs))
        if verbose
            println("a = $a, b = $b, N = $N, prior_int = $prior_integral, lh = $lh")
        end
        return lh
    catch e
        # Print out the error during optimization
        bt = catch_backtrace()
        showerror(stdout, e, bt)
        rethrow(e)
    end
end

function fit_Beta(
    X::Dict,
    n::Vector{Int64},
    ngrid::Int,
    theta0::Vector;
    lower::Vector = [0.1, 0.1],
    upper::Vector = [Inf, Inf],
    ftol::Real = 0.001
)
    """
        fit_model(X, n, ngrid, theta0, [lower, upper, ftol,])

    Maximize likelihood with respect to model parameters
    """
    # Define minimization objective
    LL(x, grad) = -loglh(x[1], x[2], X, n, ngrid)
    # Nelder Mead optimization
    opt = Opt(:LN_SBPLX, 2)
    opt.upper_bounds = upper
    opt.lower_bounds = lower
    opt.min_objective = LL
    opt.ftol_abs = ftol
    println("Optimizing....")
    (minf, minx, ret) = NLopt.optimize(opt, theta0)
    return (minf, minx, ret)
end

end
