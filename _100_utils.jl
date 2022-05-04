module Utils

import Distributions: Beta
import StatsBase: addcounts!

export loglh, loglh_truncated

function pareto_sampling(p, n)
    lambs = [[i, p[i] * n] for i in 1:length(p)];
    Q = [];
    for l in lambs
       u = rand();
       push!(Q, [l[1], (u / (1.0 - u)) / (l[2] / (1.0 - l[2]))]);
    end
    return [Int(i[1]) for i in sort(Q, by=x -> x[2])[1:n]];
end

function monte_carlo(alpha, beta, x_i, n, mcdraws)
    sum_term = 0.0;
    points = rand(Beta(alpha, beta), mcdraws);
    for p in points
        g = (n .* p).^x_i .* (1.0 .- n .* p).^(1.0 .- x_i);
        sum_term += prod(g);
    end
    return sum_term / mcdraws;
end

function lincoln(S, n)
    r = intersect(Set(S[1]), Set(S[2]));
    return prod(n) / length(r);
end

function schnabel(n, T, K)
    return sum(n[2:T] .* cumsum(n)[2:T]) / sum(values(K) .- 1.0);
end

function chao(N_o, f)
    return N_o + get(f, 1, 0)^2.0 / (2.0 * get(f, 2, 0));
end

function chao_corrected(N_o, T, f)
    m1 = 2.0 * get(f, 2, 0)  / get(f, 1, 0);
    m2 = 6.0 * get(f, 3, 0) / get(f, 1, 0);
    return (N_o
            + get(f, 1, 0)^2.0 / (2.0 * get(f, 2, 0))
            * (1.0 - m1 / T) / (1.0 - m2 / (T * m1)));
end

function jackknife(N_o, T, f, k)
    if k == 1
        return N_o + (T - 1.0) / T * get(f, 1, 0);
    elseif k == 2
        return (N_o
                + (2.0 * T - 3.0) / T * get(f, 1, 0)
                - (T - 2.0)^2.0 / (T^2.0 - T) * get(f, 2, 0));
    elseif k == 3
        return (N_o
                + (3.0 * T - 6.0) / T * get(f, 1, 0)
                - (3.0 * T^2.0 - 15.0 * T + 19.0) / (T^2.0 - T) * get(f, 2, 0)
                + (T - 3.0)^3.0 / (T * (T - 1.0) * (T - 2.0)) * get(f, 3, 0));
    elseif k == 4
        return (N_o
                + (4.0 * T - 10.0) / T * get(f, 1, 0)
                - (6.0 * T^2.0 - 36.0 * T + 55.0) / (T^2.0 - T) * get(f, 2, 0)
                + (4.0 * T^3.0 - 42.0 * T^2.0 + 148.0 * T - 175.0) / (T * (T - 1.0) * (T - 2.0)) * get(f, 3, 0)
                - (T - 4.0)^4.0 / (T * (T - 1.0) * (T - 2.0) * (T - 3.0)) * get(f, 4, 0));
    elseif k == 5
        return (N_o
                + (5.0 * T - 15.0) / T * get(f, 1, 0)
                - (10.0 * T^2.0 - 70.0 * T + 125.0) / (T^2.0 - T) * get(f, 2, 0)
                + (10.0 * T^3.0 - 120.0 * T^2.0 + 485.0 * T - 660.0) / (T * (T - 1.0) * (T - 2.0)) * get(f, 3, 0)
                - ((T - 4.0)^5.0 - (T - 5.0)^5.0) / (T * (T - 1.0) * (T - 2.0) * (T - 3.0)) * get(f, 4, 0)
                + (T - 5.0)^5.0 / (T * (T - 1.0) * (T - 2.0) * (T - 3.0) * (T - 4.0)) * get(f, 5, 0));
    else
        error("k must be between 1 and 5");
    end
end

function loglh(alpha, N_u, S, O, T, n, draws)
    N_o = length(O);
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

function loglh_truncated(alpha, N_u, S, O, T, n, draws)
    N_o = length(O);
    denom = 1 - monte_carlo(alpha, alpha * (N_u + N_o - 1), zeros(T), n, draws);
    if denom < 0
        denom = 5e-200
    end
    lh = -N_o * log(denom) + loglh(alpha, N_u, S, O, T, n, draws)
    println("alpha = $alpha, N_u = $N_u, lh = $lh");
    return lh;
end
