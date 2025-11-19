"""
Investigate the equivalence between Poisson MPLE and Horvitz-Thompson estimators
for capture-recapture data under simple random sampling.
"""

include("_100_utils.jl")

using Optim
using Distributions
using StatsBase
using Plots
using DataFrames

import .Utils

import Random: seed!

N_sim::Int64 = 10000
n_sim::Vector{Int64} = [100, 200, 300]
DRAWS::Vector{Int64} = 10:10:150
SEED::Int = 777
OUTPUT_FOLDER::String = "./_900_output/figures/appendix/"

seed!(SEED)

# Define functions for zero-truncated Poisson likelihood and MLE/MPLE estimation
function zero_trunc_pois(freqs, rate)
    out = 0.0
    for i in keys(freqs)
        out += freqs[i] * logpdf(Poisson(rate), i)
    end
    return out - sum(values(freqs)) * log(1 - pdf(Poisson(rate), 0))
end

function mle_pois(freqs)
    sum_n = sum([i * freqs[i] for i in keys(freqs)])
    res = optimize(rate -> -zero_trunc_pois(freqs, rate), 1e-8, 5 * sum_n / N_sim)
    return Optim.minimizer(res)
end

function mple_pois(freqs)
    sum_n = sum([i * freqs[i] for i in keys(freqs)])
    res = optimize(
        N -> -zero_trunc_pois(freqs, sum_n / N),
        sum(values(freqs)), 5 * N_sim
    )
    return Optim.minimizer(res)
end

# Simulate simple random sampling capture-recapture data
function simulate_capture(ns, N)
    S = []
    for n in ns
        s = sample(1:N, n, replace = false)
        push!(S, s)
    end
    K = Utils.cap_freq(S)
    return Utils.freq_of_freq(K)
end

# Function to compute mean difference between HT and MPLE estimates
function simulation_mean(ns, N, reps)
    diffs = []
    for i in 1:reps
        freqs = simulate_capture(ns, N)

        lamb_mle = mle_pois(freqs)
        N_mple = mple_pois(freqs)
        N_ht = sum(values(freqs)) / (1.0 - pdf(Poisson(lamb_mle), 0))

        delta = N_ht - N_mple
        push!(diffs, delta)
    end
    return mean(diffs)
end

# Run simulations across different sample sizes and number of capture occasions
df = DataFrame(
    n = Int[],
    t = Int[],
    delta_mean = Float64[]
)
for n in n_sim
    for t in DRAWS
        ns = repeat([n], t)
        delta_mean = simulation_mean(ns, N_sim, 1000)
        push!(df, [n, t, delta_mean])
    end
end

# Plot results
plt = plot(
    xlabel = "T",
    ylabel = "Mean difference",
    size = (1280 / 1.5, 720 / 1.5),
    leftmargin = 5Plots.mm,
    xtickfontsize = 12,
    ytickfontsize = 12,
    xguidefontsize = 12,
    yguidefontsize = 12,
    legendfontsize = 12
)

color_map = Dict(
    100 => "#003f5c",
    200 => "#bc5090",
    300 => "#ffa600"
)
for n in n_sim
    df_n = df[df.n .== n, :]
    plot!(
        plt,
        df_n.t,
        df_n.delta_mean,
        label = "n = $n",
        color = color_map[n],
        linewidth = 2
    )
end
savefig(
    plt,
    OUTPUT_FOLDER * "poisson_equivalence.pdf"
) # Figure 1 in the supplementary material of the paper
println("Plot saved to " * OUTPUT_FOLDER * "poisson_equivalence.pdf")