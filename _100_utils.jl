module Utils

import SpecialFunctions: logbeta

export simulate_caprecap


function numerator_integrand(p, alpha, N_u, N_o, S, unit)
    x_T = [unit in s_t for s_t in S];
    statistic = sum(x_T);
    y = ( p^(alpha + statistic - 1.0)
          * (1.0 - p)^(alpha * (N_o + N_u - 1.0) - 1.0)
          * prod([1.0 / length(s_t) - p for s_t in S] .^ (1 .- x_T)) );
    return y;
end

function denominator_integrand(p, alpha, N_u, N_o, S)
    y = ( p^(alpha - 1.0)
          * (1.0 - p)^(alpha * (N_o + N_u - 1.0) - 1.0)
          * prod([1.0 - length(s_t) * p for s_t in S]) );
    return y;
end

function trapezoid(f, a, b, grid)
    integral = 0.0;
    for k in 2:length(grid)
         integral += ( f(grid[k]) - f(grid[k-1]) ) * (grid[k] - grid[k-1]);
    end
    return integral * 0.5;
end
        
function likelihood(alpha, N_u, S)
    O = Set([i for j in S for i in j]);
    N_o = length(O);
    product_term = 1.0;
    for i in O
        f(x) = numerator_integrand(x, alpha, N_u, N_o, S, i);
        product_term *= trapezoid(f, 0, 1, 0.02:0.01:0.99);
    end
    return beta(alpha, alpha * (N_u + N_o - 1))^(-N_o) * product_term;
end

function loglikelihood(alpha, N_u, S)
    O = Set([i for j in S for i in j]);
    N_o = length(O);
    sum_term = 0.0;
    for i in O
        f(x) = numerator_integrand(x, alpha, N_u, N_o, S, i);
        I = trapezoid(f, 0, 1, 0.02:0.01:0.99);
        if I < 0
            I = 5e-220;
        end
        I = log(I);
        sum_term += I; 
    end
    g(x) = denominator_integrand(x, alpha, N_u, N_o, S);
    return -N_o * ( logbeta(alpha, alpha * (N_u + N_o - 1)) - trapezoid(g, 0, 1, 0.02:0.01:0.99) ) + sum_term;
end
