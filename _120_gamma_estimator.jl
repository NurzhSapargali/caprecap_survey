module GammaEstimator

using NLopt

import SpecialFunctions: digamma, loggamma

export loglh, gradient_a, gradient_Nu

function loglh(log_Nu, log_a, n, No, sum_x; verbose::Bool = true)
    Nu = exp(log_Nu)
    N = Nu + No
    alpha = exp(log_a)
    beta = N * alpha
    beta_tilde = beta + sum(n)
    sum_term = sum(loggamma.(alpha .+ sum_x))
    const_term = -No * ( log(beta_tilde^alpha - beta^alpha) + loggamma(alpha) - alpha * log(beta) )
    lh = const_term + sum_term - log(beta_tilde) * sum(sum_x)
    if verbose
        println("a = $(exp(log_a)), b = $beta, N = $(exp(log_Nu) + No), lh = $lh")
    end
    return lh
end


function gradient_a(log_Nu, log_a, n, No, sum_x)
    Nu = exp(log_Nu)
    N = Nu + No
    alpha = exp(log_a)
    beta = N * alpha
    beta_tilde = beta + sum(n)
    inner_del_num = ( beta_tilde^alpha * log(beta_tilde)
                    + beta * beta_tilde^(alpha - 1.0)
                    - beta^alpha * log(beta) )
    inner_del_denom = beta_tilde^alpha - beta^alpha
    inner_del = inner_del_num / inner_del_denom
    inner_del += digamma(alpha) - log(alpha) - 1.0
    sum_term = sum(digamma.(alpha .+ sum_x))
    out = (-No * inner_del - N / beta_tilde * sum(sum_x) + sum_term) * alpha
    println(out)
    return out
end


function gradient_Nu(log_Nu, log_a, n, No, sum_x)
    Nu = exp(log_Nu)
    N = Nu + No
    alpha = exp(log_a)
    beta = N * alpha
    beta_tilde = beta + sum(n)
    inner_del_num = alpha^2 * beta_tilde^(alpha - 1.0) - beta^(alpha - 1.0) * alpha^2
    inner_del_denom = beta_tilde^alpha - beta^alpha
    inner_del = inner_del_num / inner_del_denom
    inner_del += -alpha / N
    out = (-No * inner_del - alpha / beta_tilde * sum(sum_x)) * Nu
    println(out)
    return out
end

function fit_Gamma(theta0, n, No, sum_x; lower::Vector = [log(0.1), log(0.1)], upper::Vector = [log(10000), log(60.0)], tol = 1e-4)
    L(x) = -loglh(x[1], x[2], n, No, sum_x)
    g_Nu(x)  = -gradient_Nu(x[1], x[2], n, No, sum_x)
    g_a(x) = -gradient_a(x[1], x[2], n, No, sum_x)
#    theta = theta0
#     for i in 1:1000
#         new_Nu = theta[1] - 0.01 * g_Nu(theta)
#         new_a = theta[2] - 0.01 * g_a(theta)
#         theta = [new_Nu, new_a]
#         theta[theta .< log(0.01)] .= log(0.01)
#         println(theta)
#     end
    function objective(x::Vector, grad::Vector)
        if length(grad) > 0
            grad[1] = g_Nu(x)
            grad[2] = g_a(x)
        end
        return L(x)
    end
    opt = Opt(:LD_SLSQP, 2)
    opt.min_objective = objective
    opt.xtol_rel = tol
    opt.lower_bounds = lower
    opt.upper_bounds = upper
    opt.maxeval = 10000
    (minf,minx,ret) = optimize(opt, theta0)
#     minx = theta
#     minf = L(minx)
#     ret = 0
    return (minf, minx, ret)
end


end
