module GammaEstimator

using NLopt
using LinearAlgebra

export loglh_redux, gradient_a_redux, gradient_Nu_redux


function loglh_redux(log_alpha, log_Nu, n, ff; lamb::Real = 0.0, verbose::Bool = true)
    alpha = exp(log_alpha)
    Nu = exp(log_Nu)
    No = sum(values(ff))
    N = Nu + No
    fact_term = 0.0
    for i in keys(ff)
        multiplier = get(ff, i, 0)
        range_inc = sum([log(alpha + i - j) for j in 1:i])
        fact_term += multiplier * range_inc
    end
    ratio = 1.0 + sum(n) / alpha / N
    lh = (- sum(n) * (log(alpha) + log(N))
          - No * log(1.0 - 1.0 / (ratio^alpha))
          + fact_term
          - (No * alpha + sum(n)) * log(ratio) )
    if lamb > 0
        lh += log(lamb) - 2.0 * log(alpha) - lamb / alpha
    end
    if verbose
        println("a = $alpha, b = $(N * alpha), N = $N, lamb = $lamb, lh = $lh")
    end
    return lh
end

function gradient_prior(alpha, lamb)
    return -2.0 / alpha +  lamb / alpha^2
end

function del2_prior(alpha, lamb)
    return 2.0 / alpha^2 - 2.0 *  lamb / alpha^3
end

function inner_del_a(alpha, Nu, n, ff)
    N = Nu + sum(values(ff))
    g_a = 1.0 + sum(n) / alpha / N
    del_g_a = -sum(n) / N / alpha^2
    return (log(alpha)
            + log(N) + 1.0
            + (g_a^alpha * (log(g_a) + alpha / g_a * del_g_a) / (g_a^alpha - 1.0)) )
end

function inner_del2_a(alpha, Nu, n, ff)
    N = Nu + sum(values(ff))
    g_a = 1.0 + sum(n) / alpha / N
    del_g_a = -sum(n) / N / alpha^2
    del_g_a_a_log = (g_a^alpha * (log(g_a) + alpha / g_a * del_g_a) * log(g_a)
                     + g_a^(alpha - 1) * del_g_a)
    del_g_a_a = g_a^alpha * (log(g_a) + alpha / g_a * del_g_a)
    del_g_a_a_1 = g_a^(alpha - 1.0) * (log(g_a) + (alpha - 1.0) / g_a * del_g_a)
    del_g_a_a_a_1 = g_a^alpha + alpha * del_g_a_a - 1.0
    return (1.0 / alpha
            + (del_g_a_a_log * (g_a^alpha - 1.0) - g_a^alpha * log(g_a) * del_g_a_a)
            / (g_a^alpha - 1.0)^2
            - sum(n) / N
            * (del_g_a_a_1 * alpha * (g_a^alpha - 1.0) - g_a^(alpha - 1.0) * del_g_a_a_a_1)
            / (alpha * (g_a^alpha - 1.0))^2 )
end

function gradient_a(alpha, Nu, n, ff; lamb::Real = 0.0)
    No = sum(values(ff))
    N = Nu + No
    g_a = 1.0 + sum(n) / alpha / N
    fact_term = 0.0
    for i in keys(ff)
        multiplier = get(ff, i, 0)
        range_inc = sum([1.0 / (alpha + i - j) for j in 1:i])
        fact_term += multiplier * range_inc
    end
    gradient = (No * (log(alpha) + log(N) + 1.0 - inner_del_a(alpha, Nu, n, ff))
                + fact_term
                - sum(n) / alpha / g_a)
    if lamb > 0
        gradient += gradient_prior(alpha, lamb)
    end
    #println("....a = $alpha, b = $(N * alpha), N = $N, grad_a = $gradient")
    return gradient
end

function gradient_log_a(log_alpha, log_Nu, n, ff; lamb::Real = 0.0)
    alpha = exp(log_alpha)
    Nu = exp(log_Nu)
    return gradient_a(alpha, Nu, n, ff, lamb = lamb) * alpha
end

function del2_a(alpha, Nu, n, ff)
    No = sum(values(ff))
    N = Nu + No
    g_a = 1.0 + sum(n) / alpha / N
    fact_term = 0.0
    for i in keys(ff)
        multiplier = get(ff, i, 0)
        range_inc = sum([1.0 / (alpha + i - j)^2 for j in 1:i])
        fact_term += multiplier * range_inc
    end 
    return (No * (1.0 / alpha - inner_del2_a(alpha, Nu, n, ff))
            - fact_term
            + sum(n) / alpha / alpha / g_a^2)
end

function inner_del_N(alpha, Nu, n, ff)
    N = Nu + sum(values(ff))
    g_a = 1.0 + sum(n) / alpha / N
    return ( alpha / N - sum(n) / N / N * g_a^(alpha - 1.0) / (g_a^alpha - 1.0) )
end

function inner_del2_N(alpha, Nu, n, ff)
    N = Nu + sum(values(ff))
    g_a = 1.0 + sum(n) / alpha / N
    return ( alpha / N * (alpha - 1.0) / N * (g_a^(alpha - 2.0) - 1.0) / (g_a^alpha - 1.0)
            - (alpha / N * (g_a^(alpha - 1.0) - 1.0) / (g_a^alpha - 1.0))^2 )
end

function gradient_Nu(alpha, Nu, n, ff)
    No = sum(values(ff))
    N = Nu + No
    g_a = 1.0 + sum(n) / alpha / N
    gradient = (No * (alpha / N - inner_del_N(alpha, Nu, n, ff))
                - sum(n) / N / g_a)
    #println("....a = $alpha, b = $(N * alpha), N = $N, grad_Nu = $gradient")
    return gradient
end

function gradient_log_Nu(log_alpha, log_Nu, n, ff)
    alpha = exp(log_alpha)
    Nu = exp(log_Nu)
    return gradient_Nu(alpha, Nu, n, ff) * Nu
end

function del2_Nu(alpha, Nu, n, ff)
    No = sum(values(ff))
    N = Nu + No
    g_a = 1.0 + sum(n) / alpha / N
    return (No * (-alpha / N / N - inner_del2_N(alpha, Nu, n, ff))
            + sum(n) / N / N / g_a^2)
end

function inner_del_a_Nu(alpha, Nu, n, ff)
    N = Nu + sum(values(ff))
    g_a = 1.0 + sum(n) / alpha / N
    del_g_a = -sum(n) / N / alpha^2
    del_g_a_a_1 = g_a^(alpha - 1.0) * (log(g_a) + (alpha - 1.0) / g_a * del_g_a)
    del_g_a_a = g_a^alpha * (log(g_a) + alpha / g_a * del_g_a)
    return (1.0 / N * (g_a^(alpha - 1.0) - 1.0) / (g_a^alpha - 1.0)
            + alpha / N
            * (del_g_a_a_1 / (g_a^alpha - 1.0)
               - (g_a^(alpha - 1.0) - 1.0) * del_g_a_a / (g_a^alpha - 1.0)^2) )
end

function cross_del_a_Nu(alpha, Nu, n, ff)
    N = Nu + sum(values(ff))
    No = sum(values(ff))
    g_a = 1.0 + sum(n) / alpha / N
    del_g_a = -sum(n) / N / alpha^2
    return (No * (1.0 / N - inner_del_a_Nu(alpha, Nu, n, ff))
            + sum(n) / N / g_a^2 * del_g_a)
end

function hessian(alpha, Nu, n, ff; lamb::Real = 0.0)
    hess = zeros(2, 2)
    hess[1, 1] = del2_a(alpha, Nu, n, ff)
    if lamb > 0
        hess[1, 1] += del2_prior(alpha, lamb)
    end
    hess[2, 2] = del2_Nu(alpha, Nu, n, ff)
    hess[1, 2] = cross_del_a_Nu(alpha, Nu, n, ff)
    hess[2, 1] = hess[1, 2]
    return hess
end

function hessian_log(log_alpha, log_Nu, n, ff)
    alpha = exp(log_alpha)
    Nu = exp(log_Nu)
    H = hessian(alpha, Nu, n, ff)
    score_a = gradient_log_a(log_alpha, log_Nu, n, ff)
    score_Nu = gradient_log_Nu(log_alpha, log_Nu, n, ff)
    hess_log = zeros(2, 2)
    hess_log[1, 1] =  alpha^2 * H[1, 1] + score_a
    hess_log[2, 2] =  Nu^2 * H[2, 2] + score_Nu
    hess_log[1, 2] =  Nu * alpha * H[1, 2]
    hess_log[2, 1] =  hess_log[1, 2]
    return hess_log
end

function log_laplace(lamb, log_alpha, log_Nu, n, ff)
    alpha = exp(log_alpha)
    Nu = exp(log_Nu)
    H = hessian(alpha, Nu, n, ff, lamb = lamb)
    Hbar = det(H)
    if abs(0.5 * Hbar) < 1e-4
        Hbar  = 1e-4
    end
    ll = loglh_redux(log_alpha, log_Nu, n, ff, lamb = lamb, verbose = false) - 0.5 * log(Hbar)
    println("....a = $alpha, N = $(Nu + sum(values(ff))), lamb = $lamb, ll = $ll")
    return ll
end

function fit_Gamma(theta0, lamb, n, ff; lower::Vector = [-Inf, -Inf], upper::Vector = [10, 12], ftol = 1e-7)
    L(x) = -loglh_redux(x[1], x[2], n, ff, lamb = lamb)
    g_Nu(x)  = -gradient_log_Nu(x[1], x[2], n, ff)
    g_a(x) = -gradient_log_a(x[1], x[2], n, ff, lamb = lamb)
    function objective(x::Vector, grad::Vector)
        if length(grad) > 0
            grad[1] = g_a(x)
            grad[2] = g_Nu(x)
        end
        return L(x)
    end
    #opt = Opt(:LD_LBFGS, 2)
    opt = Opt(:LN_SBPLX, 2)
    opt.min_objective = objective
    opt.ftol_abs = ftol
    opt.lower_bounds = lower
    opt.upper_bounds = upper
    opt.maxeval = 10000
    #opt.vector_storage = 10
    (minf,minx,ret) = optimize(opt, theta0)
    return (minf, minx, ret)
end

function fit_Laplace(lamb0, log_alpha, log_Nu, n, ff)
    L(x, grad) = -log_laplace(x[1], log_alpha, log_Nu, n, ff)
    opt = Opt(:LN_SBPLX, 1)
    opt.min_objective = L
    opt.ftol_abs = 1e-7
    opt.lower_bounds = [0.15]
    opt.upper_bounds = [Inf]
    opt.maxeval = 10000
    (minf,minx,ret) = optimize(opt, [lamb0])
    return (minf, minx, ret)
end


#= function loglh_oi(log_alpha, log_Nu, log_w, n, ff, nn; verbose::Bool = true)
    w = exp(log_w)
    alpha = exp(log_alpha)
    Nu = exp(log_Nu)
    No = sum(values(ff))
    No_1 = No - sum([length(nn[i]) for i in keys(nn)])
    N = Nu + No
    b = exp(log_Nu + log(No) + log_alpha)
    log_bn = log(b + sum(n))
    coef_term = 0.0
    fact_term = 0.0
    for k in 2:maximum(keys(ff))
        multiplier = get(ff, k, 0)
        range_inc = sum([log(alpha + k - j) for j in 1:k])
        fact_term += multiplier * range_inc
        coef_term += multiplier * (alpha + k)
    end
    coef_term *= log_bn
    logratio_term = alpha * log(b) - (alpha + 1) * log_bn + log(alpha) + log(1 - w)
    last_term = 0.0
    for l in keys(nn)
        multiplier = length(get(nn, l, 0))
        last_term += multiplier * log(w + exp(logratio_term) * l)
    end
    truncation = ( log(w - 1 + (1 + sum(n) / b)^alpha)
                  - alpha * log(1 + sum(n) / b) )
    truncation *= -No
    non_one = No_1 * (log(1 - w) + alpha * log(b))
    lh = (truncation + non_one - coef_term + fact_term + last_term)
    if verbose
        println("w = $w, a = $alpha, b = $(N * alpha), N = $N, lh = $lh")
    end
    return lh
end =#
#= 
function fit_Gamma_oi(theta0, n, ff, nn; lower::Vector = [-Inf, -Inf, -Inf], upper::Vector = [10, 12, 0.0], ftol = 1e-7)
    LL(x, grad) = -loglh_oi(x[1], x[2], x[3], n, ff, nn)
    #opt = Opt(:LD_LBFGS, 2)
    opt = Opt(:LN_SBPLX, 3)
    opt.min_objective = LL
    opt.ftol_abs = ftol
    opt.lower_bounds = lower
    opt.upper_bounds = upper
    opt.maxeval = 10000
    #opt.vector_storage = 10
    (minf,minx,ret) = optimize(opt, theta0)
    return (minf, minx, ret)
end =#


end
