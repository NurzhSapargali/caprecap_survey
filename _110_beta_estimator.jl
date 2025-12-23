module BetaEstimator

import Distributions: Beta, quantile

using Optim
using StatsBase

export approx_loglike, fit_Beta

"""
    log_cond_likelihood(mu::Vector{Float64}, x::Vector{Float64})

Compute log-likelihood of a binary vector `x` given probabilities `mu`.
Uses log-space for numerical stability.
"""
function log_cond_likelihood(mu::Vector{Float64}, x::Vector{Float64})
    s = 0.0
    # Loop over capture occasions, vectorized, no bounds checking
    @inbounds @simd for i in eachindex(mu) 
        s += x[i] * log(mu[i]) + (1 - x[i]) * log1p(-mu[i])
    end
    return s
end

"""
    cond_likelihood(mu::Vector{Float64}, x::Vector{Float64})

Compute likelihood of a binary vector `x` given probabilities `mu`.
"""
function cond_likelihood(mu::Vector{Float64}, x::Vector{Float64})
    return exp(log_cond_likelihood(mu, x))
end

"""
    mc_integration(p::Vector{Float64}, x::Vector{Float64}, n::Vector{Float64})

Monte Carlo integration over capture probabilities `p` for capture history `x`
and total captures `n`.
"""
function mc_integration(p::Vector{Float64}, x::Vector{Float64}, n::Vector{Float64})
    K = length(n)
    tmp_mu = similar(n)  # preallocate once
    s = 0.0
    @inbounds for pi in p # No bounds checking
        @simd for j in 1:K # loop over capture occasions, vectorized
            tmp_mu[j] = pi * n[j]
        end
        s += cond_likelihood(tmp_mu, x)
    end
    return s / length(p)
end

"""
    approx_loglike(
        a::Real,
        Nu::Real,
        X::AbstractMatrix{<:Real},
        u::Vector{Float64};
        verbose::Bool=true
    )

Approximate log-likelihood of Beta-binomial model with parameters `a` and `Nu`
given capture history matrix `X` and Monte Carlo uniform variates `u`.

Optional arguments:
- `verbose`: if true, prints parameter values and log-likelihood (default: true)
"""
function approx_loglike(
    a::Real,
    Nu::Real,
    X::AbstractMatrix{<:Real},
    u::Vector{Float64};
    verbose::Bool = true
)
    No, K = size(X)
    N = Nu + No
    n = vec(sum(X, dims = 1))  # sample sizes per capture occasion
    n_f = Float64.(n) # convert to Float64 for arithmetic

    # Beta quantiles
    d = Beta(a, a * (N - 1.0))
    p = similar(u)
    @inbounds for i in eachindex(u) # No bounds checking
        p[i] = quantile(d, u[i])
    end
    p = p[p .<= 1.0 / maximum(n_f)]  # remove invalid quantiles

    # Compute likelihood for each row
    lh = 0.0
    @inbounds for i in 1:No # No bounds checking
        x_f = Float64.(view(X, i, :)) # convert to Float64 for arithmetic
        lh += log(mc_integration(p, x_f, n_f))
    end

    no_cap = zeros(Float64, K)
    lh += -No * log(1.0 - mc_integration(p, no_cap, n_f))

    if verbose
        println("a = $a, N = $N, lh = $lh")
    end

    return lh
end

"""
    fit_Beta(
        theta0::Vector{<:Real},
        X::AbstractMatrix{<:Real};
        lower::Vector{<:Real} = [-Inf, -Inf],
        upper::Vector{<:Real} = [Inf, Inf],
        draws::Int = 10000,
        ftol = 1e-4,
        verbose::Bool = true
    )

Fit Beta-binomial model to capture history matrix `X` (rows are individuals, columns
are capture occasions) using maximum likelihood estimation. Initial parameters
are given in `theta0` (= log(a), log(Nu)). Returns a tuple (minimum function
value, optimal parameters).

Optional arguments:
- `lower`: lower bounds for parameters (default: `[-Inf, -Inf]`)
- `upper`: upper bounds for parameters (default: `[Inf, Inf]`)
- `draws`: number of Monte Carlo draws for integration (default: 10000)
- `ftol`: function tolerance for optimization (default: 1e-4)
- `verbose`: if true, prints parameter values and log-likelihood during optimization
(default: true)
"""
function fit_Beta(
    theta0::Vector{<:Real},
    X::AbstractMatrix{<:Real},
    u::Vector{Float64}; # precomputed uniform draws for Monte Carlo integration
    lower::Vector{<:Real} = [-Inf, -Inf],
    upper::Vector{<:Real} = [Inf, Inf],
    ftol = 1e-4,
    verbose::Bool = true
)

    # Optimization function expects log-parameters (theta0 already log-transformed)
    L(x) = -approx_loglike(exp(x[1]), exp(x[2]), X, u; verbose=verbose)

    res = optimize(
        L,
        lower,
        upper,
        theta0,
        Fminbox(NelderMead()),
        Optim.Options(iterations=50, f_abstol=ftol)
    )

    return (Optim.minimum(res), Optim.minimizer(res))
end

end
