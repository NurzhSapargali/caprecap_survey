module Utils

import Distributions: Beta
import SpecialFunctions: beta

export loglh, loglh_truncated


function f(alpha_parameter, beta_parameter, 
integrand = ( prod(g) * p^(alpha_parameter - 1.0)
                      * (1.0 - p)^(beta_parameter - 1.0)
                      / beta(alpha_parameter, beta_parameter) );
function monte_carlo(alpha_parameter, beta_parameter, x_i, n, mcdraws)
    sum_term = 0.0;
    points = rand(Beta(alpha_parameter, beta_parameter), mcdraws);
    for p in points
        g = (n .* p).^x_i .* (1.0 .- n .* p).^(1.0 .- x_i);
        sum_term += prod(g);
    end
    return sum_term / mcdraws;
end

function trapezoid(alpha_parameter, beta_parameter, x_i, n, draws, a, b)
    delta = (b - a) / draws;
    points = a:delta:b;
    sum_term = 0.0;
    for p in points
        g = (n .* p).^x_i .* (1.0 .- n .* p).^(1.0 .- x_i);
        integrand = ( prod(g) * p^(alpha_parameter - 1.0)
                      * (1.0 - p)^(beta_parameter - 1.0)
                      / beta(alpha_parameter, beta_parameter) );
        if (p == a || p == b)
            sum_term += integrand;
        else
            sum_term += 2.0 * integrand;
        end
    end
    return delta / 2.0 * sum_term;
end

function loglh(alpha_parameter, N_u, S, draws)
    O = Set([i for j in S for i in j]);
    N_o = length(O);
    T = length(S);
    n = [length(i) for i in S];
    sum_term = 0.0;
    for i in O
        x_i = [i in s for s in S]; 
        I = monte_carlo(alpha_parameter, alpha_parameter * (N_o + N_u - 1), x_i, n, draws);
        #I = trapezoid(alpha_parameter, alpha_parameter * (N_o + N_u - 1), x_i, n, draws, 0.0, 1.0);
        if I < 0
            I = 5e-200;
        end
        sum_term += log(I);
    end
    #println("alpha_parameter = $alpha_parameter, N_u = $N_u, lh = $sum_term");
    return sum_term;
end

function loglh_truncated(alpha_parameter, N_u, S, draws)
    O = Set([i for j in S for i in j]);
    N_o = length(O);
    T = length(S);
    n = [length(i) for i in S];
    denom = 1 - monte_carlo(alpha_parameter, alpha_parameter * (N_u + N_o - 1), zeros(T), n, draws);
    #denom = 1 - trapezoid(alpha_parameter, alpha_parameter * (N_u + N_o - 1), zeros(T), n, draws, 0.0, 1.0);
    if denom < 0
        denom = 5e-200
    end
    lh = -N_o * log(denom) + loglh(alpha_parameter, N_u, S, draws)
    println("alpha_parameter = $alpha_parameter, N_u = $N_u, lh = $lh");
    return lh;
end
