module GammaEstimator

using NLopt

export loglh_redux, gradient_a_redux, gradient_Nu_redux


function loglh_redux(log_alpha, log_Nu, n, No, ff; verbose::Bool = true)
    alpha = exp(log_alpha)
    Nu = exp(log_Nu)
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
    if verbose
        println("a = $alpha, b = $(N * alpha), N = $N, lh = $lh")
    end
    return lh
end

function gradient_a_redux(log_alpha, log_Nu, n, No, ff)
    Nu = exp(log_Nu)
    alpha = exp(log_alpha)
    N = Nu + No
    log_aN = log_alpha + log(N)
    ratio = 1.0 + sum(n) / alpha / N
    fact_term = 0.0
    for i in keys(ff)
        multiplier = get(ff, i, 0)
        range_inc = sum([1.0 / (alpha + i - j) for j in 1:i])
        fact_term += range_inc * multiplier
    end
    inner_derivative = ( (ratio^alpha * (log_aN + log(ratio) + 1.0 / ratio) - (log_aN + 1.0) )
                        / (ratio^alpha - 1.0) )
    out = ( No * (log_aN + 1.0)
           - No * inner_derivative
           + fact_term
           - sum(n) / alpha * 1.0 / ratio )
    return out * alpha
end

function gradient_Nu_redux(log_alpha, log_Nu, n, No)
    alpha = exp(log_alpha)
    Nu = exp(log_Nu)
    N = Nu + No
    ratio = 1.0 + sum(n) / exp(log_alpha + log(N))
    inner_derivative = (ratio^(alpha - 1.0)  - 1.0) / (ratio^alpha - 1.0) * alpha / N
    out = ( alpha * No / N - No * inner_derivative - sum(n) / N * 1.0 / ratio )
    return out * Nu
end

function fit_Gamma(theta0, n, No, ff; lower::Vector = [-Inf, -Inf], upper::Vector = [10, 12], ftol = 1e-7)
    L(x) = -loglh_redux(x[1], x[2], n, No, ff)
    g_Nu(x)  = -gradient_Nu_redux(x[1], x[2], n, No)
    g_a(x) = -gradient_a_redux(x[1], x[2], n, No, ff)
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
