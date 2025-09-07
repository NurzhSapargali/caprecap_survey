"""
Utility functions for capture–recapture simulations and sampling.
"""
module Utils

using StatsBase
using Distributions

export ar_pareto_sample, write_row, freq_of_freq,
       cap_freq, read_captures, lower_pr, simulate_samples,
       create_folder_if_not_exists

"""
    write_row(filename::String, row::Vector)

Append `row` as a comma-separated line to `filename`.
If the file does not exist, it will be created.
"""
function write_row(filename::String, row::Vector)
    open(filename, "a") do io
        for (i,j) in enumerate(row)
            if i != length(row)
                print(io, j, ",")
            else
                print(io, j, "\n")
            end
        end
    end
end

"""
    create_folder_if_not_exists(folder::String)

Create the folder if it does not already exist.
"""
function create_folder_if_not_exists(folder::String)
    if !isdir(folder)
        mkpath(folder)
    end
end

function read_captures(filename::String)
    raw = open(f -> read(f, String), filename);
    out = Vector{Int64}[];
    for l in split(raw, '\n')
        if length(l) != 0
            push!(out, parse.(Int, split(l, ",")));
        end
    end
    return out;
end

"""
    pareto_sample(p::Vector{Float64}, n::Int64)

Draw `n` indices from a population with size probabilities `p` using Pareto sampling.
"""
function pareto_sample(p::Vector{Float64}, n::Int64)
    pr = n * p # Inclusion probabilities
    N = length(pr)
    U = rand(N)
    Q = (U ./ (1.0 .- U)) ./ (pr ./ (1.0 .- pr)) # Transformed uniform variables
    return sortperm(Q)[1:n] # Return indices of the smallest n Q values
end

"""
    ar_pareto_sample(p::Vector{Float64}, n::Int64)

Draw `n` indices from a population with unequal inclusion probabilities `p`
using the acceptance–rejection Pareto sampling algorithm (Bondesson et. al, 2006).

References
----------
Bondesson, L., Traat, I. and Lundqvist, A. (2006), Pareto Sampling versus
Sampford and Conditional Poisson Sampling. Scandinavian Journal of Statistics,
33: 699-720. https://doi.org/10.1111/j.1467-9469.2006.00497.x
"""
function ar_pareto_sample(p::Vector{Float64}, n::Int64)
    pr = n * p
    N = length(pr)
    d = sum(pr .* (1.0 .- pr)) # page 701
    Jk = 1.0 .+ pr .* (pr .- 0.5) / d # Approximation 1 on page 702
    A = minimum(Jk) # Algorithm 1 on page 709

    S = []
    accept = false
    while (!accept)
        S = pareto_sample(p, n)
        car = sum(1.0 .- pr[S]) / sum(Jk[S] / A .* (1.0 .- pr[S])) # Equation (22) on page 710
        U = rand()
        #println("car = $car, U = $U")
        if (U <= car)
            accept = true
        end
    end
    return S
end

"""
    simulate_samples(pop::Int64, T::Int64, alpha::Float64, r::Float64, q::Float64)

Simulate capture–recapture samples for a population of size `pop` over `T` occasions
with heterogeneity parameter `alpha` and Negative Binomial sample size parameters `r` and `
q`.
"""
function simulate_samples(pop::Int64, T::Int64, alpha::Float64, r::Float64, q::Float64, T_min::Int64 = 5)
    n = rand(Distributions.NegativeBinomial(r, q), T) .+ 1 # Sample sizes for T occasions

    # Generate proportional to size probabilities from truncated Beta distribution
    b = alpha * (pop - 1.0)
    d = Distributions.truncated(Distributions.Beta(alpha, b), upper = 1.0 / maximum(n))
    p = rand(d, pop)

    # Generate samples using adaptive rejection Pareto sampling
    sum_n_T_min = 0
    samples = []
    O = Set{Int64}()

    # Ensure that there is non-zero recaptures in the first T_min samples
    while !(length(O) < sum_n_T_min)
        samples = [ar_pareto_sample(p, i) for i in n]
        O = Set([i for j in samples[1:T_min] for i in j])
        sum_n_T_min = sum([length(s) for s in samples[1:T_min]])
    end

    return samples
end

"""
    cap_freq(S::Vector)

Compute the capture frequencies of individuals in the list of samples `S`.
Returns a dictionary mapping individual IDs to their capture counts.
"""
function cap_freq(S::Vector)
    K = Dict{Int, Int}();
    for s in S
        addcounts!(K, s);
    end
    return K;
end

"""
    freq_of_freq(K::Dict)

Compute the frequency of frequencies from the capture frequency dictionary `K`.
Returns a dictionary mapping capture counts to the number of individuals with that count.
"""
function freq_of_freq(K::Dict)
    return countmap(values(K));
end

"""
    resample_f(f::Dict, N_hat, B::Int)

Resample the frequency of frequencies `f` using a multinomial distribution
with estimated population size `N_hat` for `B` bootstrap samples.
Returns a vector of `B` resampled frequency of frequencies dictionaries.
"""
function resample_f(f::Dict, N_hat, B::Int)
    No = sum(values(f))
    vals_f = collect(values(f))
    
    p_hat = vals_f / N_hat
    push!(p_hat, (N_hat - No) / N_hat) # Probability of unobserved individuals

    out = Vector{Dict{Int, Int}}(undef, B)
    for b in 1:B
        res = rand(Distributions.Multinomial(Int(round(N_hat)), p_hat))
        fb = Dict(
            collect(keys(f))[i] => res[i] for i in 1:length(vals_f) if res[i] > 0
        )
        out[b] = fb
    end
    return out
end


function lower_pr(N, n, q)
    return exp((n - 1) * ( log(N * (1.0 - q) - 1.0) - log(N) - log(1.0 - q) ) + log(n - 1.0 + N * (1.0 - q)) - log(N) - log(1.0 - q))
end

end
