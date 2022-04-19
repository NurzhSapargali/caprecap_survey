module Utils

import Distributions: Beta

export loglh, loglh_truncated


function monte_carlo(alpha, beta, x_i, n, mcdraws)
    sum_term = 0.0;
    points = rand(Beta(alpha, beta), mcdraws);
    for p in points
        g = ((n .* p).^(-1.0) .- 1).^(1 .- x_i);
        sum_term += prod(g);
    end
    return sum_term / mcdraws;
end

function loglh(alpha, N_u, S, mcdraws)
    O = Set([i for j in S for i in j]);
    N_o = length(O);
    T = length(S);
    n = [length(i) for i in S];
    first_sum_term = 0.0;
    second_sum_term = 0.0;
    for t in 1:T
        first_sum_term += log(alpha + T - t) - log(alpha * (N_o + N_u) + T - t);
    end
    for i in O
        x_i = [i in s for s in S]; 
        I = monte_carlo(alpha + T, alpha * (N_o + N_u - 1), x_i, n, mcdraws);
        second_sum_term += log(I);
    end
    return N_o * first_sum_term + second_sum_term;
end

function loglh_truncated(alpha, N_u, S, mcdraws)
    O = Set([i for j in S for i in j]);
    N_o = length(O);
    T = length(S);
    n = [length(i) for i in S];
    product_term = 1.0;
    sum_term = 0.0;
    for t in 1:T
        product_term *= (alpha * (N_u + N_o) + T - t) / (alpha + T - t);
    end
    product_term *= prod(n);
    product_term -= monte_carlo(alpha + T, alpha * (N_u + N_o - 1), zeros(T), n, mcdraws);
    for i in O
        x_i = [i in s for s in S];
        I = monte_carlo(alpha + T, alpha * (N_o + N_u - 1), x_i, n, mcdraws);
        sum_term += log(I);
    end
    return N_o * log(product_term) + sum_term;
end
