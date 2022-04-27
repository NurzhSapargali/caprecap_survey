module Utils

import Distributions: Beta
import StatsBase: addcounts!

export loglh, loglh_truncated

function pareto_sampling(p, n)
    lambs = [[i, p[i] * n] for i in 1:length(p)];
    Q = [[i[1], i[2] / sum(lambs)[2]] for i in lambs];
    for i in 1:N
       u = rand();
       Q[i][2] = (u / (1.0 - u)) / (Q[i][2] / (1.0 - Q[i][2]));
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

function lincoln(S)
    r = intersect(Set(S[1]), Set(S[2]));
    return prod([length(s) for s in S]) / length(r);
end

function schnabel(S)
    K = Dict{Float64, Int}();
    for s in S
        addcounts!(K, s);
    end
    n = [length(s) for s in S];
    T = length(S);
    return sum(n[2:T] .* cumsum(n)[2:T]) / sum(values(K) .- 1.0);
end

function chao(S)
    K = Dict{Float64, Int}();
    for s in S
        addcounts!(K, s);
    end
    f = countmap(values(K));
    O = Set([i for j in S for i in j]);
    return length(O) + f[1]^2.0 / (2.0 * f[2]);
end

function chao_corrected(S)
    K = Dict{Float64, Int}();
    for s in S
        addcounts!(K, s);
    end
    f = countmap(values(K));
    O = Set([i for j in S for i in j]);
    m1 = 2.0 * f[2] / f[1];
    m2 = 6.0 * f[3] / f[1];
    return length(O) + f[1]^2.0 / (2.0 * f[2]) * (1.0 - m1 / length(S)) / (1.0 - m2 / (length(S) * m1));
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
