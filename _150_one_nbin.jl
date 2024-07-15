module OneNbin

using Distributions
using StatsBase
using SpecialFunctions
using Optim

export log_nbin_trunc, log_likelihood, fit_oi_nbin_trunc, fit_oi_geom_trunc


function log_nbin_trunc(y, a, N, sum_n)
    ratio = 1.0 + sum_n / a / N
    return ( -loggamma(y + 1) - log(ratio^a - 1.0)
             + y * log(sum_n) - y * log(ratio)
             - y * log(a) - y * log(N)
             + sum([log(y + a - i) for i in 1:y])
            )
end

function log_likelihood(logit_w, log_a, log_Nu, f; verbose = true)
    No = sum(values(f))
    N = exp(log_Nu) + No
    a = exp(log_a)
    w = 1.0 / (1.0 + exp(-logit_w))
    sum_n = sum([get(f, i, 0) * i for i in keys(f)])
    singles = f[1]
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

function gradient_w(logit_w, log_a, log_Nu, f)
    No = sum(values(f))
    N = exp(log_Nu) + No
    a = exp(log_a)
    w = 1.0 / (1.0 + exp(-logit_w))
    sum_n = sum([get(f, i, 0) * i for i in keys(f)])
    singles = f[1]
    single_prob = exp(log_nbin_trunc(1, a, N, sum_n))
    score_w = ( singles / (w + (1.0 - w) * single_prob)
                * (1.0 - single_prob)
                - (No - singles) / (1.0 - w)
               )
    return score_w * w * (1.0 - w)
end


function log_nbin_gradient_a(y, a, N, sum_n)
    ratio = 1.0 + sum_n / a / N
    return ( -ratio^a / (ratio^a - 1.0) * (log(ratio) - (ratio - 1.0) / ratio)
             + y / ratio * (ratio - 1.0) / a
             - y / a
             + sum([1.0 / (y + a - i) for i in 1:y]) 
             )
end


function gradient_log_a(logit_w, log_a, log_Nu, f)
    No = sum(values(f))
    N = exp(log_Nu) + No
    a = exp(log_a)
    w = 1.0 / (1.0 + exp(-logit_w))
    sum_n = sum([get(f, i, 0) * i for i in keys(f)])
    singles = f[1]
    single_prob = exp(log_nbin_trunc(1, a, N, sum_n))
    score_a = ( singles / (w + (1.0 - w) * single_prob)
                * (1.0 - w) * log_nbin_gradient_a(1, a, N, sum_n) * single_prob
                + sum([f[i] * log_nbin_gradient_a(i, a, N, sum_n) for i in keys(f) if i > 1])
               )
    return score_a * a
end


function log_nbin_gradient_Nu(y, a, N, sum_n)
    ratio = 1.0 + sum_n / a / N
    return ( ratio^(a - 1.0) / (ratio^a - 1.0) * sum_n / N / N
             + y / ratio * (ratio - 1.0) / N
             - y / N )
end


function gradient_log_Nu(logit_w, log_a, log_Nu, f)
    No = sum(values(f))
    N = exp(log_Nu) + No
    a = exp(log_a)
    w = 1.0 / (1.0 + exp(-logit_w))
    sum_n = sum([get(f, i, 0) * i for i in keys(f)])
    singles = f[1]
    single_prob = exp(log_nbin_trunc(1, a, N, sum_n))
    score_Nu = ( singles / (w + (1.0 - w) * single_prob)
                 * (1.0 - w) * log_nbin_gradient_Nu(1, a, N, sum_n) * single_prob
                 + sum([f[i] * log_nbin_gradient_Nu(i, a, N, sum_n) for i in keys(f) if i > 1])
                )
    return score_Nu * exp(log_Nu)
end

function gradient_lamb_r(lamb_r, d)
    a = d / (d - 1.0)^2
    return (1.0 / lamb_r + 1.0 / (1.0 - lamb_r) - a / (1.0 - lamb_r)^2)
end

function fit_oi_nbin_trunc(
    theta, f;
    ftol = 1e-5, lower = [-Inf, -Inf, -Inf], upper = [Inf, 10.0, 23.0]
)
    L(x) = -log_likelihood(x[1], x[2], x[3], f)
    function g!(G::Vector, x::Vector)
        G[1] = -gradient_w(x[1], x[2], x[3], f)
        G[2] = -gradient_log_a(x[1], x[2], x[3], f)
        G[3] = -gradient_log_Nu(x[1], x[2], x[3], f)
    end
    res = optimize(L, g!, lower, upper, theta, Fminbox(GradientDescent()), Optim.Options(iterations = 50, f_tol = ftol))
    (minf, minx) = (Optim.minimum(res), Optim.minimizer(res))
    return (minf, minx)
end


function fit_oi_geom_trunc(
    theta, f;
    ftol = 1e-5, lower = [-Inf, -Inf], upper = [Inf, 23.0]
)
    L(x) = -log_likelihood(x[1], log(1.0), x[2], f)
    function g!(G::Vector, x::Vector)
        G[1] = -gradient_w(x[1], log(1.0), x[2], f)
        G[2] = -gradient_log_Nu(x[1], log(1.0), x[2], f)
    end
    res = optimize(L, g!, lower, upper, theta, Fminbox(GradientDescent()), Optim.Options(iterations = 50, f_tol = ftol))
    (minf, minx) = (Optim.minimum(res), Optim.minimizer(res))
    return (minf, minx)
end

end