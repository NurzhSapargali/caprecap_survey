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
seed!(777)

OUTPUT_FOLDER::String = "./_900_output/figures/appendix/"

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

function simulate_capture(ns, N)
    S = []
    for n in ns
        s = sample(1:N, n, replace = false)
        push!(S, s)
    end
    K = Utils.cap_freq(S)
    return Utils.freq_of_freq(K)
end

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

plt = plot(
    xlabel = "T",
    ylabel = "Mean difference",
    size = (1280 / 1.5, 720 / 1.5),
    leftmargin = 5Plots.mm
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
        color = color_map[n]
    )
end
savefig(
    plt,
    OUTPUT_FOLDER * "poisson_equivalence.pdf"
)