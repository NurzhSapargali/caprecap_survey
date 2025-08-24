include("../_100_utils.jl")

using .Utils

using Test

using Distributions
using StatsBase
using HypothesisTests

import Random: seed!

seed!(42)
ALPHAS::Vector{Float64} = [0.5, 1.0, 2.0]
N::Int64 = 1000
n::Int64 = 30
REPS::Int64 = 50000

println("KS tests between true inclusion probabilities and relative frequencies from AR-Pareto sampling")
println("---------------------------------------------------------------")
for a in ALPHAS
    println("Testing alpha = $a, N = $N, n = $n, REPS = $REPS")

    b = a * (N - 1.0)
    d = Distributions.truncated(Distributions.Beta(a, b), upper = 1.0 / n)
    p = rand(d, N)

    samples = [Utils.ar_pareto_sample(p, n) for _ in 1:REPS]
    freq = Utils.cap_freq(samples)
    rel_freq = collect(values(freq)) ./ REPS
    pr = n * p
    ks = HypothesisTests.ApproximateTwoSampleKSTest(pr, rel_freq)

    @test pvalue(ks) > 0.75
    println("p-value = ", pvalue(ks))
    println("-----")
end