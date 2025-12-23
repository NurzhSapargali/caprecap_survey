"""
OneBin module for fitting zero-truncated one-inflated Negative Binomial models.
"""
module OneNbin

using Distributions
using StatsBase
using SpecialFunctions
using Optim

import LogExpFunctions: logistic, logit

export fit_oi_nbin_trunc, fit_oi_geom_trunc, w_hat

"""
    log_nbin(y, a, N, sum_n)

Log-probability of observing count `y` from a Negative Binomial distribution
with heterogeneity parameter `a`, population size `N`, and total sample size `sum_n`.
"""
function log_nbin(y, a, N, sum_n)
    return ( -SpecialFunctions.loggamma(y + 1)
             + a * (log(a) + log(N) - log(sum_n + a * N))
             + y * (log(sum_n) - log(sum_n + a * N))
             + SpecialFunctions.loggamma(y + a) - SpecialFunctions.loggamma(a) )
end

"""
    log_nbin_trunc(y, a, N, sum_n)

Log-probability of observing count `y` from a zero-truncated Negative Binomial
distribution with heterogeneity parameter `a`, population size `N`, and total
sample size `sum_n`.
"""
function log_nbin_trunc(y, a, N, sum_n)
    ratio = 1.0 + sum_n / (a * N)
    return ( -SpecialFunctions.loggamma(y + 1) - log(ratio^a - 1.0)
             + y * (log(sum_n) - log(sum_n + a * N)) 
             + SpecialFunctions.loggamma(y + a) - SpecialFunctions.loggamma(a) )
end

"""
    log_likelihood(logit_w, log_a, log_Nu, f; verbose = true)

Log-likelihood of the one-inflated zero-truncated Negative Binomial model
with parameters `logit_w` (logit of one-inflation probability), `log_a`
(log of heterogeneity parameter), and `log_Nu` (log of the number of unobserved
individuals), given the frequency of frequencies `f`. If `verbose` is true,
prints the parameter values and log-likelihood.
"""
function log_likelihood(logit_w, log_a, log_Nu, f; verbose = true)
    No = sum(values(f)) # Number of observed individuals across all capture occasions
    N = exp(log_Nu) + No
    a = exp(log_a)
    w = logistic(logit_w)

    sum_n = sum([get(f, i, 0) * i for i in keys(f)])
    singles = get(f, 1, 0)
    single_prob = exp(log_nbin_trunc(1, a, N, sum_n))

    lh = ( singles * log(w + (1.0 - w) * single_prob)
           + (No - singles) * log(1.0 - w)
           + sum([f[i] * log_nbin_trunc(i, a, N, sum_n) for i in keys(f) if i > 1])
        )
    if verbose
        println("w = $w, a = $a, N = $N, lh = $lh")
    end
    return lh
end

"""
    log_likelihood_ztoi(w, log_a, log_Nu, f; verbose = true)

Log-likelihood of the zero-truncated one-inflated Negative Binomial model
with parameters `w` (one-inflation probability), `log_a` (log of
heterogeneity parameter), and `log_Nu` (log of the number of unobserved individuals),
given the frequency of frequencies `f`.
If `verbose` is true, prints the parameter values and log-likelihood.
"""
function log_likelihood_ztoi(w, log_a, log_Nu, f; verbose = true)
    No = sum(values(f)) # Number of observed individuals across all capture occasions
    N = exp(log_Nu) + No
    a = exp(log_a)

    sum_n = sum([get(f, i, 0) * i for i in keys(f)])
    singles = get(f, 1, 0)
    single_prob = exp(log_nbin(1, a, N, sum_n))
    zero_prob = exp(log_nbin(0, a, N, sum_n))

    lh_single = (
        -log(1.0 - (1.0 - w) * zero_prob)
        + log(1.0 - (1.0 - w) * (1.0 - single_prob))
    )
    lh_nonsingles = sum(
        [
            f[i] * (log(1.0 - w)
            - log(1.0 - (1.0 - w) * zero_prob)
            + log_nbin(i, a, N, sum_n)) for i in keys(f) if i > 1
        ]
    )

    lh = singles * lh_single + lh_nonsingles
    if verbose
        println("w = $w, log_a = $log_a, log_Nu = $log_Nu, N = $N, lh = $lh")
    end
    return lh
end  

"""
    w_hat(log_a, log_Nu, f)

Compute the maximum likelihood estimate of the one-inflation probability `w`
given parameters `log_a` (log of heterogeneity parameter) and `log_Nu`
and frequency of frequencies `f`.
"""

function w_hat(log_a, log_Nu, f)
    No = sum(values(f)) # Number of observed individuals across all capture occasions
    N = exp(log_Nu) + No
    a = exp(log_a)

    sum_n = sum([get(f, i, 0) * i for i in keys(f)])
    singles = get(f, 1, 0)
    single_prob = exp(log_nbin_trunc(1, a, N, sum_n))

    w_hat = (singles - single_prob * No) / (No - No * single_prob)
    return clamp(w_hat, 1e-6, 1.0 - 1e-6)
end

"""
    profile_log_likelihood(log_a, log_Nu, f; verbose = true)

Profile log-likelihood of the zero-truncated one-inflated Negative Binomial model
with parameters `log_a` (log of heterogeneity parameter) and `log_Nu` (log of
the number of unobserved individuals), given the frequency of frequencies `f`.
The one-inflation probability `w` is profiled out. If `verbose` is true, prints
the parameter values and log-likelihood.
"""
function profile_log_likelihood(log_a, log_Nu, f; verbose = true)
    No = sum(values(f)) # Number of observed individuals across all capture occasions
    N = exp(log_Nu) + No
    a = exp(log_a)

    w = w_hat(log_a, log_Nu, f)

    lh = log_likelihood(logit(w), log_a, log_Nu, f; verbose = false)
    if verbose
        println("w = $w, log_a = $log_a, log_Nu = $log_Nu, N = $N, lh = $lh")
    end
    return lh
end

"""
    gradient_w(logit_w, log_a, log_Nu, f)

Gradient of the log-likelihood with respect to `logit_w`.
"""
function gradient_w(logit_w, log_a, log_Nu, f)
    No = sum(values(f))
    N = exp(log_Nu) + No
    a = exp(log_a)
    w = logistic(logit_w)

    sum_n = sum([get(f, i, 0) * i for i in keys(f)])
    singles = get(f, 1, 0)
    single_prob = exp(log_nbin_trunc(1, a, N, sum_n))

    score_w = ( singles / (w + (1.0 - w) * single_prob)
                * (1.0 - single_prob)
                - (No - singles) / (1.0 - w) )
    return score_w * w * (1.0 - w)
end


"""
    ratio_derivative(a, N, sum_n)

Derivative of the ratio term used in gradients with respect to `a`.
"""
function ratio_derivative(a, N, sum_n)
    ratio = 1.0 + sum_n / a / N
    return ratio^a * (log(ratio) - sum_n / a / N * ratio^(-1.0))
end

"""
    log_nbin_gradient_a(y, a, N, sum_n)

Gradient of the log-probability of observing count `y` from a zero-truncated
Negative Binomial distribution with respect to the heterogeneity parameter `a`.
"""
function log_nbin_gradient_a(y, a, N, sum_n)
    ratio = 1.0 + sum_n / a / N
    return ( -y * N  / (sum_n + a * N)
             - ratio_derivative(a, N, sum_n) / (ratio^a - 1.0)
            + SpecialFunctions.digamma(y + a) - SpecialFunctions.digamma(a) )
end

"""
    gradient_log_a(log_a, log_Nu, f)

Gradient of the profile log_likelihood with respect to log of `a`
"""
function gradient_log_a(log_a, log_Nu, f)
    No = sum(values(f))
    N = exp(log_Nu) + No
    a = exp(log_a)
    
    sum_n = sum([get(f, i, 0) * i for i in keys(f)])
    singles = get(f, 1, 0)
    single_prob = exp(log_nbin_trunc(1, a, N, sum_n))

    w_hat = (singles - single_prob * No) / (No - No * single_prob)
    w_hat = clamp(w_hat, 1e-6, 1.0 - 1e-6)

    score_a = ( singles / (w_hat + (1.0 - w_hat) * single_prob)
                * (1.0 - w_hat) * log_nbin_gradient_a(1, a, N, sum_n) * single_prob
                + sum([f[i] * log_nbin_gradient_a(i, a, N, sum_n) for i in keys(f) if i > 1]) )
    return score_a * a
end

"""
    log_nbin_gradient_Nu(y, a, N, sum_n)

Gradient of the log-probability of observing count `y` from a zero-truncated
Negative Binomial distribution with respect to the number of unobserved individuals
`Nu`.
"""
function log_nbin_gradient_Nu(y, a, N, sum_n)
    ratio = 1.0 + sum_n / a / N
    return ( -y * a / (sum_n + a * N)
             + ratio^(a - 1.0) / (ratio^a - 1.0) * (sum_n / N^2) )
end

"""
    gradient_log_Nu(log_a, log_Nu, f)

Gradient of the profile log_likelihood with respect to log of `Nu`
"""
function gradient_log_Nu(log_a, log_Nu, f)
    No = sum(values(f))
    N = exp(log_Nu) + No
    a = exp(log_a)

    sum_n = sum([get(f, i, 0) * i for i in keys(f)])
    singles = get(f, 1, 0)
    single_prob = exp(log_nbin_trunc(1, a, N, sum_n))

    w_hat = (singles - single_prob * No) / (No - No * single_prob)
    w_hat = clamp(w_hat, 1e-6, 1.0 - 1e-6)

    score_Nu = ( singles / (w_hat + (1.0 - w_hat) * single_prob)
                 * (1.0 - w_hat) * log_nbin_gradient_Nu(1, a, N, sum_n) * single_prob
                 + sum([f[i] * log_nbin_gradient_Nu(i, a, N, sum_n) for i in keys(f) if i > 1]) )
    return score_Nu * exp(log_Nu)
end

"""
    fit_oi_nbin_trunc(
        theta, f;
        frel_tol = 1e-5,
        lower = [-Inf, -750.0],
        upper = [10.0, 23.0],
        verbose = true,
        method = Optim.LBFGS()
    )

Fit the zero-truncated one-inflated Negative Binomial model to frequency of
frequencies `f` using initial parameter estimates `theta` (a vector of `[log_a,
log_Nu]`).

Optional arguments:
- `f_reltol`: function tolerance for optimization (default: 1e-5)
- `lower`: lower bounds for parameters (default: `[-Inf, -750.0]`)
- `upper`: upper bounds for parameters (default: `[10.0, 23.0]`)
- `verbose`: if true, prints parameter values and log-likelihood during optimization
(default: true)
- `method`: optimization method from Optim.jl (default: Optim.LBFGS())

Returns a tuple `(minf, minx)` where `minf` is the minimum negative log-likelihood
and `minx` is the vector of estimated parameters.
"""
function fit_oi_nbin_trunc(
    theta, f;
    f_reltol = 1e-5,
    lower = [-Inf, -750.0],
    upper = [10.0, 23.0],
    verbose = true,
    method = Optim.LBFGS()
)
    # Objective function
    L(x) = -profile_log_likelihood(x[1], x[2], f; verbose = verbose)
    
    # Gradients
    function g!(G::Vector, x::Vector)
        G[1] = -gradient_log_a(x[1], x[2], f)
        G[2] = -gradient_log_Nu(x[1], x[2], f)
    end

    # Optimize using box-constrained method
    res = Optim.optimize(
        L,
        g!,
        lower,
        upper,
        theta,
        Optim.Fminbox(method),
        Optim.Options(f_reltol = f_reltol)
    )
    (minf, minx) = (Optim.minimum(res), Optim.minimizer(res))

    return (minf, minx)
end

"""
    fit_oi_geom_trunc(
        theta, f;
        f_reltol = 1e-5,
        lower = [-Inf],
        upper = [23.0],
        verbose = true,
        method = Optim.LBFGS()
    )

Fit the zero-truncated one-inflated Geometric model to frequency of frequencies `f`
using initial parameter estimate `theta` (a vector of `[log_Nu]`).

Optional arguments:
- `f_reltol`: function tolerance for optimization (default: 1e-5)
- `lower`: lower bounds for parameters (default: `[-Inf]`)
- `upper`: upper bounds for parameters (default: `[23.0]`)
- `verbose`: if true, prints parameter values and log-likelihood during optimization
(default: true)
- `method`: optimization method from Optim.jl (default: Optim.LBFGS())

Returns a tuple `(minf, minx)` where `minf` is the minimum negative log-likelihood
and `minx` is the vector of estimated parameters.
"""
function fit_oi_geom_trunc(
    theta, f;
    f_reltol = 1e-5,
    lower = [-750.0],
    upper = [23.0],
    verbose = true,
    method = Optim.LBFGS()
)
    # Objective function
    L(x) = -profile_log_likelihood(log(1.0), x[1], f; verbose = verbose)
    
    # Gradients
    function g!(G::Vector, x::Vector)
        G[1] = -gradient_log_Nu(log(1.0), x[1], f)
    end

    # Optimize using box-constrained method
    res = Optim.optimize(
        L,
        g!,
        lower,
        upper,
        theta,
        Optim.Fminbox(method),
        Optim.Options(f_reltol = f_reltol)
    )
    (minf, minx) = (Optim.minimum(res), Optim.minimizer(res))

    return (minf, minx)
end


"""
    fit_trunc_nbin_oi(
        theta, f;
        freltol = 1e-5,
        lower = [0.0, -Inf, -750.0],
        upper = [1.0, 10.0, 23.0],
        verbose = true,
        method = Optim.NelderMead()
    )

Fit the zero-truncated Negative Binomial one-inflated model to frequency of frequencies `f`
using initial parameter estimates `theta` (a vector of `[w, log_a, log_Nu]`).
Optional arguments:
- `f_reltol`: function tolerance for optimization (default: 1e-5
- `lower`: lower bounds for parameters (default: `[0.0, -Inf, -750.0]`)
- `upper`: upper bounds for parameters (default: `[1.0, 10.0, 23.0]`)
- `verbose`: if true, prints parameter values and log-likelihood during optimization
(default: true)
- `method`: optimization method from Optim.jl (default: Optim.NelderMead())
Returns a tuple `(minf, minx)` where `minf` is the minimum negative log-likelihood
and `minx` is the vector of estimated parameters.
"""
function fit_trunc_nbin_oi(
    theta, f;
    freltol = 1e-5,
    lower = [1e-5, -Inf, -750.0],
    upper = [0.999, 10.0, 23.0],
    verbose = true,
    method = Optim.NelderMead()
)
    # Objective function
    L(x) = -log_likelihood_ztoi(x[1], x[2], x[3], f; verbose = verbose)

    # Optimize using box-constrained method
    res = Optim.optimize(
        L,
        lower,
        upper,
        theta,
        Optim.Fminbox(method),
        Optim.Options(f_reltol = freltol)
    )
    (minf, minx) = (Optim.minimum(res), Optim.minimizer(res))

    return (minf, minx)
end

end