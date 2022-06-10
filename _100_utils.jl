module Utils

import Distributions: Beta
using StatsBase

export sampford_sample, lincoln, schnabel,
       chao, chao_corrected, jackknife,
       loglh_truncated, write_row
       
function write_row(filename, row)
    open(filename, "a") do io
        for (i,j) in enumerate(row)
            if i != length(row)
                print(io, j, ",");
            else
                print(io, j, "\n");
            end
        end
    end
end

function sampford_sample(p, n)
    pr = n * p;
    pr = Dict(enumerate(pr));
    s = Set();
    while length(s) != n
        key = sample(1:length(p), Weights(p));
        val = pop!(pr, key);
        new_pr = collect(values(pr));
        new_pr = new_pr ./ (1.0 .- new_pr);
        new_pr = new_pr ./ sum(new_pr);
        s = Set(sample(collect(keys(pr)), Weights(new_pr), n));
        pr[key] = val;
    end
    println("Generated")
    return collect(s);
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

function schnabel(S, n)
    m = 0;
    N = 0;
    pool = Set();
    denom = [];
    num = [];
    R = Set();
    for t in 2:length(S)
        pool = union(pool, S[t-1]);
        u = n[t-1] - length(R);
        R = intersect(S[t], pool);
        m += u;
        push!(num, n[t] * m);
        push!(denom, length(R));
    end
    return sum(num) / sum(denom);
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

end

