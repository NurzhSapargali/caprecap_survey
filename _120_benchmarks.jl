module Benchmarks

import SpecialFunctions: beta

using LinearAlgebra
using NLopt

export lincoln, schnabel, chao,
       zelterman, jackknife, huggins,
       turing, conway_maxwell, turing_geometric


function lincoln(S::Vector,
                 n::Vector{Int64})
    r = intersect(Set(S[1]), Set(S[2]));
    return prod(n) / length(r);
end

function schnabel(S::Vector,
                  n::Vector{Int64})
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

function chao(N_o::Int64, f::Dict{Int64, Int64})
    return N_o + get(f, 1, 0)^2.0 / (2.0 * get(f, 2, 0));
end

function zelterman(N_o::Int64, f::Dict{Int64, Int64})
    return N_o / (1.0 - exp(-2.0 * get(f, 2, 0) / get(f, 1, 0)));
end

function jackknife(N_o::Int64, T::Int64, f::Dict{Int64, Int64}, k::Int64)
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

function conway_maxwell(N_o::Int64, f::Dict{Int64, Int64})
    y = [log(i * get(f, i, 0) / get(f, i - 1, 0)) for i in 2:length(f)];
    if ((Inf in y) || (-Inf in y))
        return "Nonsequential frequencies";
    end
    x = [log(i) for i in 2:length(f)];
    x = hcat(ones(length(x)), x);
    w = diagm([1.0 / get(f, i - 1, 0) + 1.0 / get(f, i, 0) for i in 2:length(f)]);
    w = inv(w);
    if det(transpose(x) * w * x) == 0
        return "Noninvertible X'WX";
    end
    b = inv(transpose(x) * w * x) * transpose(x) * w * y;
    N = N_o + get(f, 1, 0) * exp(-b[1]);
    return N;
end

function moments_huggins(theta::Vector{Float64},
                         p::Vector{Float64},
                         T::Int64)
    a = theta[1];
    b = theta[2];
    cdf = 1.0 -  beta(a, b + T) / beta(a, b);
    first = a / (a + b) /  cdf;
    second = (a * b * (a + b + T) / (T * (a + b + T) * (a + b)^2) + (a / (a + b))^2) / cdf;
    out = [first, second];
    return out .- p;
end

function huggins(T::Int64,
                 K::Dict{Int64, Int64})
    param = [sum(values(K)) / length(K) / T, sum(values(K) .^ 2) / length(K) / (T ^ 2)];
    f(x, grad) = sum(moments_huggins(x, param, T)) ^ 2;
    opt = Opt(:LN_SBPLX, 2);
    opt.lower_bounds = [0.01, 0.01];
    opt.min_objective = f;
    (minf, minx, ret) = NLopt.optimize(opt, [5.0, 5.0]);
    mu = (values(K) .+ minx[1]) / (sum(minx) + T);
    return sum(1.0 ./ (1.0 .- (1.0 .- mu).^T));
end

function turing_geometric(N_o::Int64, f::Dict{Int64, Int64}, T::Int64)
    denom = sum([i * get(f, i, 0) for i in 1:maximum(keys(f))]);
    return N_o / (1.0 - sqrt(get(f, 1, 0) / denom));
end

function turing(N_o::Int64, f::Dict{Int64, Int64}, T::Int64)
    denom = sum([i * get(f, i, 0) for i in 1:maximum(keys(f))]);
    return N_o / (1.0 - (get(f, 1, 0) / denom) ^ (T / (T - 1)));
end

end
