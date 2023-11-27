module Benchmarks

import SpecialFunctions: beta

using LinearAlgebra
using NLopt
using RCall

export lincoln, schnabel, chao,
       zelterman, jackknife, morgan_ridout,
       turing, conway_maxwell, turing_geometric


function lincoln(S::Vector,
                 n::Vector{Int})
    r = intersect(Set(S[1]), Set(S[2]));
    return prod(n) / length(r);
end

function schnabel(S::Vector,
                  n::Vector{Int})
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

function chao(N_o::Int64, f::Dict)
    return N_o + get(f, 1, 0)^2.0 / (2.0 * get(f, 2, 0));
end

function zelterman(N_o::Int64, f::Dict)
    return N_o / (1.0 - exp(-2.0 * get(f, 2, 0) / get(f, 1, 0)));
end

function jackknife(N_o::Int64, T::Int64, f::Dict, k::Int64)
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

function conway_maxwell(N_o::Int64, f::Dict)
    y = [log(i * get(f, i, 0) / get(f, i - 1, 0)) for i in 2:maximum(keys(f))]
    x = [log(i) for i in 2:maximum(keys(f))]
    x = hcat(ones(length(x)), x)
    faults = findall((y .== -Inf) .| (y .== Inf) .| (isnan.(y)))
    inds = setdiff(1:length(y), faults)
    y = y[inds,:]
    x = x[inds,:]
    w = diagm([1.0 / get(f, i - 1, 0) + 1.0 / get(f, i, 0) for i in 2:maximum(keys(f))])
    w = w[inds, inds]
    w = inv(w)
    if length(y) < 2
        return -999
    else
        b = inv(transpose(x) * w * x) * transpose(x) * w * y
        N = N_o + get(f, 1, 0) * exp(-b[1])
        return N
    end
end

function chao_lee_jeng(N_o::Int64, f::Dict, T::Int, n::Vector, bias::Int = 0)
    seq = 1:T
    cdenom = 0
    gamma_num = 0
    double_nsum = 0
    for k in 1:T
        cdenom += k * get(f, k, 0)
        gamma_num += k * (k - 1) * get(f, k, 0)
        remain = seq[seq .> k]
        for j in remain
            double_nsum += n[k] * n[j]
        end
    end
    if bias == 1
        C = 1.0 - (get(f, 1, 0) - 2.0 * get(f, 2, 0) / (T - 1)) / cdenom
    elseif bias == 2
        C = 1.0 - (get(f, 1, 0) - 2.0 * get(f, 2, 0) / (T - 1) + 6.0 * get(f, 3, 0) / ((T - 1) * (T - 2))) / cdenom
    else
        C = 1.0 - get(f, 1, 0) / cdenom
    end
    No_hat = N_o / C 
    gamma = maximum([No_hat * gamma_num / (2.0 * double_nsum) - 1.0, 0.0])^2
    return No_hat + get(f, 1, 0) / C * gamma
end
        
    
function morgan_ridout(f::Dict, T::Int, filename::String; seed::Int = 1)
    skinks = [get(f, i, 0) for i in 1:maximum(keys(f))]
    hat = -999
    R"""
    set.seed($seed)
    source($filename)
    out <- 0.0
    fail <- TRUE
    tol <- 10
    counter <- 0
    while (fail){
        tryCatch(
            out <- capture.output(fitonemodel($skinks, $T, model="binbbin")),
            fail <<- FALSE,
            error = function(e){ fail <<- TRUE }
        )
        if (counter > tol){
            fail <<- FALSE
        }
        counter <- counter + 1
    }
    """
    out = @rget out
    if (length(out) == 1)
        hat = -999
    else
        out = out[occursin.("Estimated total population", out)][1]
        out = strip(split(out, "= ")[2])
        hat = parse(Float64, out)
    end
    return hat
end

function turing_geometric(N_o::Int, f::Dict, T::Int)
    denom = sum([i * get(f, i, 0) for i in 1:maximum(keys(f))]);
    return N_o / (1.0 - sqrt(get(f, 1, 0) / denom));
end

function turing(N_o::Int, f::Dict, T::Int)
    denom = sum([i * get(f, i, 0) for i in 1:maximum(keys(f))]);
    return N_o / (1.0 - (get(f, 1, 0) / denom) ^ (T / (T - 1)));
end

end
