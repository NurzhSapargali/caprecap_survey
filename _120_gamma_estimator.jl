module GammaEstimator

using NLopt
using LinearAlgebra

export loglh
export gradient_a, gradient_log_a, gradient_Nu, gradient_log_Nu
export hessian, hessian_log
export fit_Gamma


function loglh(log_alpha, log_Nu, n, ff; verbose::Bool = true)
    alpha = exp(log_alpha)
    Nu = exp(log_Nu)
    No = sum(values(ff))
    N = Nu + No
    fact_term = 0.0
    ratio = 1.0 + sum(n) / alpha / N
    # upper_N = 2 * No - get(ff, maximum(keys(ff)), 0)
    for i in keys(ff)
        multiplier = get(ff, i, 0)
        range_inc = sum([log(alpha + i - j) for j in 1:i])
        fact_term += multiplier * range_inc
    end
    lh = ( -No * log(ratio^alpha - 1.0)
           - sum(n) * (log(alpha) + log(N) + log(ratio))
           + fact_term )
    if verbose
        println("a = $alpha, b = $(N * alpha), N = $N, lh = $lh")
    end
    return lh
end


function gradient_a(alpha, Nu, n, ff)
    No = sum(values(ff))
    N = Nu + No
    ratio = 1.0 + sum(n) / alpha / N
    ratio_derivative = -sum(n) / alpha^2 / N
    fact_term = 0.0
    for i in keys(ff)
        multiplier = get(ff, i, 0)
        range_inc = sum([1.0 / (alpha + i - j) for j in 1:i])
        fact_term += multiplier * range_inc
    end
    gradient = ( -No / (ratio^alpha - 1.0)
                 * (ratio^alpha * (log(ratio) + alpha / ratio * ratio_derivative))
                 - sum(n) * (1.0 / alpha + 1.0 / ratio * ratio_derivative)
                 + fact_term )
    #println("....a = $alpha, b = $(N * alpha), N = $N, grad_a = $gradient")
    return gradient
end


function gradient_log_a(log_alpha, log_Nu, n, ff)
    alpha = exp(log_alpha)
    Nu = exp(log_Nu)
    return gradient_a(alpha, Nu, n, ff) * alpha
end


function gradient_Nu(alpha, Nu, n, ff)
    No = sum(values(ff))
    N = Nu + No
    ratio = 1.0 + sum(n) / alpha / N
    ratio_derivative = -sum(n) / alpha / N^2
    gradient = ( -No / (ratio^alpha - 1.0) * alpha * ratio^(alpha - 1.0) * ratio_derivative
                 - sum(n) * (1.0 / N + ratio_derivative / ratio) )
    #println("....a = $alpha, b = $(N * alpha), N = $N, grad_Nu = $gradient")
    return gradient
end


function gradient_log_Nu(log_alpha, log_Nu, n, ff)
    alpha = exp(log_alpha)
    Nu = exp(log_Nu)
    return gradient_Nu(alpha, Nu, n, ff) * Nu
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


function inner_del2_N(alpha, Nu, n, ff)
    N = Nu + sum(values(ff))
    g_a = 1.0 + sum(n) / alpha / N
    return ( alpha / N * (alpha - 1.0) / N * (g_a^(alpha - 2.0) - 1.0) / (g_a^alpha - 1.0)
            - (alpha / N * (g_a^(alpha - 1.0) - 1.0) / (g_a^alpha - 1.0))^2 )
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


function fit_Gamma(theta0,
    n,
    ff;
    lower::Vector = [-Inf, -Inf],
    upper::Vector = [10, 12],
    ftol = 1e-7
)
    L(x) = -loglh(x[1], x[2], n, ff)
    g_Nu(x)  = -gradient_log_Nu(x[1], x[2], n, ff)
    g_a(x) = -gradient_log_a(x[1], x[2], n, ff)
    function objective(x::Vector, grad::Vector)
        if length(grad) > 0
            grad[1] = g_a(x)
            grad[2] = g_Nu(x)
        end
        return L(x)
    end
    opt = Opt(:LD_LBFGS, 2)
    #opt = Opt(:LN_SBPLX, 2)
    opt.min_objective = objective
    opt.ftol_abs = ftol
    opt.lower_bounds = lower
    opt.upper_bounds = upper
    opt.maxeval = 10000
    #opt.vector_storage = 10
    (minf,minx,ret) = optimize(opt, theta0)
    return (minf, minx, ret)
end

end