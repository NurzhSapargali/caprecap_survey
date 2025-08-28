module Benchmarks

using LinearAlgebra
using RCall

export lincoln, schnabel, chao,
       zelterman, jackknife, morgan_ridout,
       turing, conway_maxwell, turing_geometric


"""    
    lincoln(S::Vector, n::Vector{Int})

Lincoln-Petersen estimator for population size given two capture samples `S` and their sizes `n`.
"""
function lincoln(S::Vector, n::Vector{Int})
    r = intersect(Set(S[1]), Set(S[2]))
    return prod(n) / length(r)
end

"""
    schnabel(S::Vector, n::Vector{Int})

Schnabel estimator for population size given capture samples `S` and their sizes `n`.
"""
function schnabel(S::Vector, n::Vector{Int})
    pool = Set(S[1]) # Set of marked individuals prior to capture t
    num = 0.0
    denom = 0.0
    for t in 2:length(S)
        m = length([s for s in S[t] if s in pool]) # Number of recaptures in sample t
        denom += m
        num += length(pool) * n[t] 
        union!(pool, Set(S[t])) # Update pool with individuals from sample t
    end
    return num / denom
end

"""
    chao(f::Dict)

Chao estimator for population size given frequency of capture frequencies `f`.
"""
function chao(f::Dict)
    return sum(values(f)) + get(f, 1, 0)^2.0 / (2.0 * get(f, 2, 0))
end

"""
    zelterman(f::Dict)

Zelterman estimator for population size given frequency of capture frequencies `f`.
"""
function zelterman(f::Dict)
    return sum(values(f)) / (1.0 - exp(-2.0 * get(f, 2, 0) / get(f, 1, 0)))
end

"""
    jackknife(T::Int, f::Dict, k::Int)

Jackknife estimator of order `k` for population size given number of capture occasions `T`
and frequency of capture frequencies `f`.
"""
function jackknife(T::Int, f::Dict, k::Int)
    N_o = sum(values(f))
    if k == 1
        return N_o + (T - 1.0) / T * get(f, 1, 0)
    elseif k == 2
        return (N_o
                + (2.0 * T - 3.0) / T * get(f, 1, 0)
                - (T - 2.0)^2.0 / (T^2.0 - T) * get(f, 2, 0))
    elseif k == 3
        return (N_o
                + (3.0 * T - 6.0) / T * get(f, 1, 0)
                - (3.0 * T^2.0 - 15.0 * T + 19.0) / (T^2.0 - T) * get(f, 2, 0)
                + (T - 3.0)^3.0 / (T * (T - 1.0) * (T - 2.0)) * get(f, 3, 0))
    elseif k == 4
        return (N_o
                + (4.0 * T - 10.0) / T * get(f, 1, 0)
                - (6.0 * T^2.0 - 36.0 * T + 55.0) / (T^2.0 - T) * get(f, 2, 0)
                + (4.0 * T^3.0 - 42.0 * T^2.0 + 148.0 * T - 175.0) / (T * (T - 1.0) * (T - 2.0)) * get(f, 3, 0)
                - (T - 4.0)^4.0 / (T * (T - 1.0) * (T - 2.0) * (T - 3.0)) * get(f, 4, 0))
    elseif k == 5
        return (N_o
                + (5.0 * T - 15.0) / T * get(f, 1, 0)
                - (10.0 * T^2.0 - 70.0 * T + 125.0) / (T^2.0 - T) * get(f, 2, 0)
                + (10.0 * T^3.0 - 120.0 * T^2.0 + 485.0 * T - 660.0) / (T * (T - 1.0) * (T - 2.0)) * get(f, 3, 0)
                - ((T - 4.0)^5.0 - (T - 5.0)^5.0) / (T * (T - 1.0) * (T - 2.0) * (T - 3.0)) * get(f, 4, 0)
                + (T - 5.0)^5.0 / (T * (T - 1.0) * (T - 2.0) * (T - 3.0) * (T - 4.0)) * get(f, 5, 0))
    else
        error("k must be between 1 and 5")
    end
end

"""
    conway_maxwell(f::Dict)

Ratio regression estimator for Conway-Maxwell-Poisson distributed counts 
given frequency of capture frequencies `f`.
"""
function conway_maxwell(f::Dict)
    N_o = sum(values(f))
    # Prepare response vector and design matrix
    y = [log(i * get(f, i, 0) / get(f, i - 1, 0)) for i in 2:maximum(keys(f))]
    x = log.(2:maximum(keys(f)))
    x = hcat(ones(length(x)), x)

    # Remove infinite and NaN values from y and corresponding rows from x
    faults = findall((y .== -Inf) .| (y .== Inf) .| (isnan.(y)))
    inds = setdiff(1:length(y), faults)
    y = y[inds,:]
    x = x[inds,:]

    # Weight matrix for weighted least squares
    w = LinearAlgebra.inv(
        LinearAlgebra.diagm(
            [1.0 / get(f, i - 1, 0) + 1.0 / get(f, i, 0) for i in 2:maximum(keys(f))]
        )
    )
    w = w[inds, inds] # Adjust weight matrix to match filtered y and x

    # If there are not enough data points to perform regression, return -999
    if length(y) < 2
        return -999
    else
        b = inv(transpose(x) * w * x) * transpose(x) * w * y
        N = N_o + get(f, 1, 0) * exp(-b[1])
        return N
    end
end


"""
    chao_lee_jeng(f::Dict, T::Int, n::Vector, bias::Int = 0)

Chao-Lee-Jeng estimator for population size given frequency of capture frequencies
`f`, number of capture occasions `T`, and sample sizes `n`.
If `bias` is 1 or 2, applies bias correction of the first or second order, respectively.
"""
function chao_lee_jeng(f::Dict, T::Int, n::Vector, bias::Int = 0)
    N_o = sum(values(f))
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
    if bias == 0
        C = 1.0 - get(f, 1, 0) / cdenom
    elseif bias == 1
        C = 1.0 - (get(f, 1, 0) - 2.0 * get(f, 2, 0) / (T - 1)) / cdenom
    elseif bias == 2
        C = 1.0 - (get(f, 1, 0) - 2.0 * get(f, 2, 0) / (T - 1) + 6.0 * get(f, 3, 0) / ((T - 1) * (T - 2))) / cdenom
    else
        error("bias must be 0, 1, or 2")
    end
    No_hat = N_o / C 
    gamma2 = maximum([No_hat * gamma_num / (2.0 * double_nsum) - 1.0, 0.0])
    return No_hat + get(f, 1, 0) / C * gamma2
end

"""
    morgan_ridout(f::Dict, T::Int, filename::String; seed::Int = 1)

Binomial and beta-binomial mixture estimator for population size given frequency
of capture frequencies `f` and number of capture occasions `T`. Requires R script
implemented by M.S.Ridout at `filename` to perform estimation.
The optional `seed` parameter sets the random seed in R for reproducibility.
"""
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
    # Retry fitting the model if it fails, up to `tol` times
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

"""
    turing_geometric(f::Dict)

Turing estimator for population size under assumption of geometric counts given frequency
of capture frequencies `f`
"""
function turing_geometric(f::Dict)
    denom = sum([i * get(f, i, 0) for i in 1:maximum(keys(f))]);
    return sum(values(f)) / (1.0 - sqrt(get(f, 1, 0) / denom));
end

"""
    turing(f::Dict, T::Int)

Turing estimator for population size of homogeneous population given frequency
of capture frequencies `f` and number of capture occasions `T`.
"""
function turing(f::Dict, T::Int)
    denom = sum([i * get(f, i, 0) for i in 1:maximum(keys(f))]);
    return sum(values(f)) / ( 1.0 - (get(f, 1, 0) / denom)^(T / (T - 1)) );
end

end
