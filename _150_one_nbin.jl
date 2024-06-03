using Distributions
using StatsBase
using SpecialFunctions
using NLopt


function log_prior_d(d, lamb_r)
    lamb = lamb_r / (1.0 - lamb_r)
    dist = Exponential(1.0 / lamb)
    return logpdf(dist, d) - 2.0 * log(d)
end

function log_nbin_trunc(y, a, N, sum_n)
    ratio = 1.0 + sum_n / a / N
    return ( -loggamma(y + 1) - log(ratio^a - 1.0)
             + y * log(sum_n) - y * log(ratio)
             - y * log(a) - y * log(N)
             + sum([log(y + a - i) for i in 1:y])
            )
end

function log_likelihood(w, log_a, log_Nu, f; verbose = true)
    No = sum(values(f))
    N = exp(log_Nu) + No
    a = exp(log_a)
    #w = 1.0 / (1.0 + exp(-logit_w))
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

function log_joint(f, d, w, log_Nu, lamb_r; verbose = true)
    a = d / (d - 1.0)^2
    log_a = log(a)
    log_prior = log_prior_d(d, lamb_r)
    log_lh = log_likelihood(w, log_a, log_Nu, f; verbose = verbose)
    return log_prior + log_lh
end

function gradient_w(w, log_a, log_Nu, f)
    No = sum(values(f))
    N = exp(log_Nu) + No
    a = exp(log_a)
    #w = 1.0 / (1.0 + exp(-logit_w))
    sum_n = sum([get(f, i, 0) * i for i in keys(f)])
    singles = f[1]
    single_prob = exp(log_nbin_trunc(1, a, N, sum_n))
    score_w = ( singles / (w + (1.0 - w) * single_prob)
                * (1.0 - single_prob)
                - (No - singles) / (1.0 - w)
               )
    return score_w #* w * (1.0 - w)
end


function log_nbin_gradient_a(y, a, N, sum_n)
    ratio = 1.0 + sum_n / a / N
    return ( -ratio^a / (ratio^a - 1.0) * (log(ratio) - (ratio - 1.0) / ratio)
             + y / ratio * (ratio - 1.0) / a
             - y / a
             + sum([1.0 / (y + a - i) for i in 1:y]) 
             )
end


function gradient_log_a(w, log_a, log_Nu, f)
    No = sum(values(f))
    N = exp(log_Nu) + No
    a = exp(log_a)
    #w = 1.0 / (1.0 + exp(-logit_w))
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


function gradient_log_Nu(w, log_a, log_Nu, f)
    No = sum(values(f))
    N = exp(log_Nu) + No
    a = exp(log_a)
    #w = 1.0 / (1.0 + exp(-logit_w))
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

function grid_posterior(f, w, log_Nu, lamb_r; grid_delta = 0.01)
    grid_vals = grid_delta:grid_delta:(1.0 - grid_delta)
    log_p = [log_joint(f, i, w, log_Nu, lamb_r; verbose = false) for i in grid_vals]
    c = 700.0 - maximum(log_p)
    app_p = exp.(log_p .+ c)
    return Dict(grid_vals[i] => (app_p / sum(app_p))[i] for i in eachindex(grid_vals))
end

function sample_grid(f, w, log_Nu, lamb_r, n; grid_delta = 0.01)
    grid_post = grid_posterior(f, w, log_Nu, lamb_r; grid_delta = grid_delta)
    p = collect(values(grid_post))
    return sample(collect(keys(grid_post)), pweights(p), n)
end

function fit_oi_nbin_trunc(
    theta, f;
    ftol = 1e-6, lower = [0, -Inf, -Inf], upper = [1.0, 9.9, 23.0]
)
    function objective(x::Vector, grad::Vector)
        if length(grad) > 0
            grad[1] = -gradient_w(x[1], x[2], x[3], f)
            grad[2] = -gradient_log_a(x[1], x[2], x[3], f)
            grad[3] = -gradient_log_Nu(x[1], x[2], x[3], f)
        end
        return -log_likelihood(x[1], x[2], x[3], f)
    end
    opt = Opt(:LD_LBFGS, 3)
    opt.lower_bounds = lower
    opt.upper_bounds = upper
    opt.min_objective = objective
    opt.ftol_rel = ftol
    (minf, minx, ret) = optimize(opt, theta)
    return (minf, minx, ret)
end


function fit_oi_geom_trunc(
    theta, f;
    ftol = 1e-6, lower = [0.0, -Inf], upper = [1.0, 23.0]
)
    function objective(x::Vector, grad::Vector)
        if length(grad) > 0
            grad[1] = -gradient_w(x[1], log(1.0), x[2], f)
            grad[2] = -gradient_log_Nu(x[1], log(1.0), x[2], f)
        end
        return -log_likelihood(x[1], log(1.0), x[2], f)
    end
    opt = Opt(:LD_LBFGS, 2)
    opt.lower_bounds = lower
    opt.upper_bounds = upper
    opt.min_objective = objective
    opt.ftol_rel = ftol
    (minf, minx, ret) = optimize(opt, theta)
    return (minf, minx, ret)
end

function m_step(
    f, w, log_Nu, lamb_r, post_d;
    xtol = 1e-3, lower = [0, -Inf, 0.001], upper = [1.0, 23.0, 0.999]
)
    post_a = post_d ./ (1.0 .- post_d)
    function objective(x::Vector, grad::Vector)
        if length(grad) > 0
            grad[1] = -mean([gradient_w(x[1], log(i), x[2], f) for i in post_a])
            grad[2] = -mean([gradient_log_Nu(x[1], log(i), x[2], f) for i in post_a])
            grad[3] = -mean([gradient_lamb_r(x[3], i) for i in post_d])
        end
        out = -mean([log_joint(f, i, x[1], x[2], x[3]; verbose = false) for i in post_d])
        println("w = $(x[1]), N = $(exp(x[2]) + sum(values(f))), lamb_r = $(x[3]), out = $(-out)")
        return out
    end
    theta = [w, log_Nu, lamb_r]
    opt = Opt(:LN_SBPLX, 3)
    opt.lower_bounds = lower
    opt.upper_bounds = upper
    opt.min_objective = objective
    opt.xtol_rel = xtol
    (minf, minx, ret) = optimize(opt, theta)
    return (minf, minx, ret)
end
old_lh = 0.0
counter = 0.0
converged = false
while !converged
    post_d = sample_grid(f, w, log_Nu, lamb_r, 5000; grid_delta = 0.0001);
    (minf, minx, ret) = m_step(f, w, log_Nu, lamb_r, post_d)
    (lh, (w, log_Nu, lamb_r), status) = m_step(f, w, log_Nu, lamb_r, post_d)
    lh_delta = abs(lh - old_lh)
    if lh_delta < 1e-4 && counter > 10
        converged = true
    end
    old_lh = lh
    counter += 1
    println("Iteration $counter")
end