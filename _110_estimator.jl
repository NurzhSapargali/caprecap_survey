module Estimator

using NLopt

import SpecialFunctions: digamma, loggamma

export loglh, gradient_a, gradient_Nu, update_a, update_Nu

function loglh(log_Nu, log_a, n, No, sum_x)
    Nu = exp(log_Nu)
    N = Nu + No
    alpha = exp(log_a)
    beta = N * alpha
    beta_tilde = beta + sum(n)
    sum_term = sum(loggamma.(alpha .+ sum_x))
    const_term = -No * ( log(beta_tilde^alpha - beta^alpha) + loggamma(alpha) - alpha * log(beta) )
    return const_term + sum_term - log(beta_tilde) * sum(sum_x)
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
    return (-No * inner_del - N / beta_tilde * sum(sum_x) + sum_term) * alpha
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
    return (-No * inner_del - alpha / beta_tilde * sum(sum_x)) * Nu
end

# function update_a(log_a, log_Nu, n, No, sum_x; eta = 0.01)
#     return log_a + eta * gradient_a(log_Nu, log_a, n, No, sum_x)
# end
# 
# function update_Nu(log_Nu, log_a, n, No, sum_x; eta = 0.5)
#     return log_Nu + eta * gradient_Nu(log_Nu, log_a, n, No, sum_x)
# end
# 
# function coordinate_ascent(theta0, n, No, sum_x; eta = 0.1, tol = 0.001, maxiter = 10000)
#     log_ahat = theta0[1]
#     log_Nuhat = theta0[2]
#     converged = false
#     iters = 0
#     while !converged
#         iters += 1
#         old_Nuhat = log_Nuhat
#         old_ahat = log_ahat
#         log_Nuhat = Estimator.update_Nu(log_Nuhat, log_ahat, n, No, sum_x; eta = 0.1)
#         log_ahat = Estimator.update_a(log_ahat, log_Nuhat, n, No, sum_x; eta = 0.1)
#         if sum(isnan.([log_ahat, log_Nuhat])) > 0.0
#             return [old_ahat, old_Nuhat]
#         end
#         if ((abs(log_Nuhat - old_Nuhat) < tol)&(abs(log_ahat - old_ahat) < tol))|(maxiter <= iters)
#             converged = true
#         end
#         println(iters)
#     end
#     return [log_ahat, log_Nuhat]
# end

function fit_model(theta0, n, No, sum_x; tol = 1e-4)
    L(x) = -loglh(x[1], x[2], n, No, sum_x)
    g_Nu(x)  = -gradient_Nu(x[1], x[2], n, No, sum_x)
    g_a(x) = -gradient_a(x[1], x[2], n, No, sum_x)
    function objective(x::Vector, grad::Vector)
        if length(grad) > 0
            grad[1] = g_Nu(x)
            grad[2] = g_a(x)
        end
        #println(x)
        return L(x)
    end
    opt = Opt(:LD_MMA, 2)
    opt.min_objective = objective
    opt.xtol_rel = tol
    opt.maxeval = 10000
    (minf,minx,ret) = optimize(opt, theta0)
    return minx
end


end
