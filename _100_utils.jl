module Utils

import Distributions: Beta
import StatsBase: addcounts!

export loglh, loglh_truncated


function monte_carlo(alpha, beta, x_i, n, mcdraws)
    sum_term = 0.0;
    points = rand(Beta(alpha, beta), mcdraws);
    for p in points
        g = (n .* p).^x_i .* (1.0 .- n .* p).^(1.0 .- x_i);
        sum_term += prod(g);
    end
    return sum_term / mcdraws;
end

function lincoln(S)
    r = intersect(Set(S[1]), Set(S[2]));
    return prod([length(s) for s in S]) / length(r);
end

function schnabel(S)
    K = Dict{Int, Int}();
    for s in S
        addcounts!(K, s);
    end
    n = [length(s) for s in S];
    T = length(S);
    return sum(n[2:T] .* cumsum(n)[2:T]) / sum(values(K) .- 1);
end

function loglh(alpha, N_u, S, draws)
    O = Set([i for j in S for i in j]);
    N_o = length(O);
    T = length(S);
    n = [length(i) for i in S];
    sum_term = 0.0;
    for i in O
        x_i = [i in s for s in S]; 
        I = monte_carlo(alpha, alpha * (N_o + N_u - 1), x_i, n, draws);
        if I < 0
            I = 5e-200;
        end
        sum_term += log(I);
    end
    return sum_term;
end

function loglh_truncated(alpha, N_u, S, draws)
    O = Set([i for j in S for i in j]);
    N_o = length(O);
    T = length(S);
    n = [length(i) for i in S];
    denom = 1 - monte_carlo(alpha, alpha * (N_u + N_o - 1), zeros(T), n, draws);
    if denom < 0
        denom = 5e-200
    end
    lh = -N_o * log(denom) + loglh(alpha, N_u, S, draws)
    println("alpha = $alpha, N_u = $N_u, lh = $lh");
    return lh;
end
