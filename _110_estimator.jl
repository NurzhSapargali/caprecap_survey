module Estimator

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
    fraction_num = ( beta_tilde^alpha * (log(beta_tilde) + beta / beta_tilde)
                    - beta^alpha * (log(beta) + 1.0) )
    fraction_denom = beta_tilde^alpha - beta^alpha
    fraction = fraction_num / fraction_denom
    fraction += digamma(alpha) - log(beta) - 1.0
    sum_term = sum(digamma.(alpha .+ sum_x))
    return (-No * fraction - N / beta_tilde * sum(sum_x) + sum_term)
end


function gradient_Nu(log_Nu, log_a, n, No, sum_x)
    Nu = exp(log_Nu)
    N = Nu + No
    alpha = exp(log_a)
    beta = N * alpha
    beta_tilde = beta + sum(n)
    fraction_num = alpha^2.0 * ( beta_tilde^(alpha - 1.0) - beta^(alpha - 1.0) )
    fraction_denom = beta_tilde^alpha - beta^alpha
    fraction = fraction_num / fraction_denom
    fraction += -alpha / N
    return (-No * fraction - alpha / beta_tilde * sum(sum_x))
end

function update_a(a, log_Nu, n, No, sum_x; eta = 0.01)
    return a + eta * gradient_a(log_Nu, log(a), n, No, sum_x)
end

function update_Nu(Nu, log_a, n, No, sum_x; eta = 0.5)
    return Nu + eta * gradient_Nu(log(Nu), log_a, n, No, sum_x)
end

end
