module GammaEstimator

using NLopt

export loglh, gradient_a, gradient_Nu

function loglh(Nu, alpha, n, No, x_sums; verbose::Bool = true)
    N = Nu + No
    fact_term = 0.0
    for i in keys(x_sums)
        incs = get(x_sums, i, 0)
        range_inc = sum([log(alpha + incs - j) for j in 1:incs])
        fact_term += range_inc
    end
    large_term = alpha * log(alpha * N + sum(n)) + log( 1.0 - ( alpha * N / (alpha * N + sum(n)) )^alpha )
    lh = ( No * alpha * (log(alpha) + log(N))
          - log(alpha * N + sum(n)) * sum(n)
          + fact_term
          - No * large_term )
    if verbose
        println("a = $alpha, b = $(N * alpha), N = $N, lh = $lh")
    end
    return lh
end
    
function loglh_redux(Nu, alpha, n, No, x_sums; verbose::Bool = true)
    N = Nu + No
    fact_term = 0.0
    for i in keys(x_sums)
        incs = get(x_sums, i, 0)
        range_inc = sum([log(alpha + incs - j) for j in 1:incs])
        fact_term += range_inc
    end
    ratio = sum(n) / alpha / N
    lh = (- sum(n) * (log(alpha) + log(N))
          - No * log(1.0 - 1.0 / ((1.0 + ratio)^alpha))
          + fact_term
          - (No * alpha + sum(n)) * log(1.0 + ratio) )
    if verbose
        println("a = $alpha, b = $(N * alpha), N = $N, lh = $lh")
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

function fit_Gamma(theta0, n, No, x_sums; lower::Vector = [0.1, 0.1], upper::Vector = [Inf, Inf], xtol = 1e-7)
    L(x, grad) = -loglh_redux(x[1], x[2], n, No, x_sums)
    opt = Opt(:LN_SBPLX, 2)
    opt.min_objective = L
    opt.xtol_rel = xtol
    opt.lower_bounds = lower
    opt.upper_bounds = upper
    opt.maxeval = 10000
    (minf,minx,ret) = optimize(opt, theta0)
    return (minf, minx, ret)
end


end
