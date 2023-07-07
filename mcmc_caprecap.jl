import SpecialFunctions: loggamma
import Distributions: Beta, truncated, Gamma, Geometric, Poisson, Normal, Bernoulli, logpdf

using StatsBase
using Plots
using Turing

function loglike(sum_x, sum_n, No, a, Nu)
    a_tilde = sum_x + a
    N = No + Nu
    b = a * N + sum_n
    return (loggamma(a_tilde)
            - loggamma(a)
            + a * log(a)
            + a * log(N)
            - a_tilde * log(b))
end

function loga_prior(a, c, d)
    d = Gamma(c, 1.0 / d)
    return logpdf(d, a)
end

function logNu_prior(Nu, sigma2)
    d = Normal(0, sigma2)
    return logpdf(d, Nu) - log(Nu)
end

function log_posterior(a, Nu, X, sum_n, c, d, sigma2)
    unadj = sum([loglike(X[i], sum_n, length(X), a, Nu) for i in keys(X)])
    renorm = -length(X) * log(1.0 - exp(loglike(0, sum_n, length(X), a, Nu)))
    comp_loglh = unadj + renorm
    return comp_loglh + loga_prior(a, c, d) + logNu_prior(Nu, sigma2)
end

function mcmc(a0, Nu0, X, sum_n, c, d, p; chain = 1000, prop_var_Nu = 1.0, prop_var_a = 1.0)
    old_a = a0
    old_Nu = Nu0
    Nu_trace = [old_Nu]
    a_trace = [old_a]
    accept_count = 0
    for t in 1:chain
        new_u = rand(Normal(log(old_a), prop_var_a))
        new_a = exp(new_u)
        new_z = rand(Normal(log(old_Nu), prop_var_Nu))
        new_Nu = exp(new_z)
        old_post = log_posterior(old_a, old_Nu, X, sum_n, c, d, p)
        new_post = log_posterior(new_a, new_Nu, X, sum_n, c, d, p)
        old_prop = logpdf(Normal(new_u, prop_var_a), log(old_a)) - log(old_a) + logpdf(Normal(new_z, prop_var_Nu), log(old_Nu)) - log(old_Nu)
        new_prop = logpdf(Normal(log(old_a), prop_var_a), new_u) - new_u + logpdf(Normal(log(old_Nu), prop_var_Nu), new_z) - new_z
        numerator = new_post + old_prop
        denom = old_post + new_prop
        accept = minimum([1, exp(numerator - denom)])
        if (rand(Bernoulli(accept)) == 1)
            old_a = new_a
            old_Nu = new_Nu
            accept_count += 1
        end
        push!(Nu_trace, old_Nu)
        push!(a_trace, old_a)
    end
    println(accept_count / chain)
    return hcat(a_trace, Nu_trace)
end

function summarize(theta, No, burn_in)
    out1 = [quantile(theta[burn_in:size(theta)[1], 1], 0.025),
           median(theta[burn_in:size(theta)[1], 1]),
           quantile(theta[burn_in:size(theta)[1], 1], 0.975)]
    out2 = [quantile(theta[burn_in:size(theta)[1], 2], 0.025),
            median(theta[burn_in:size(theta)[1], 2]),
            quantile(theta[burn_in:size(theta)[1], 2], 0.975)]
    return hcat(out1, out2 .+ No)
end

function plot_results(theta, No, burn_in)
    p1 = plot(theta[:,1], label = "alpha")
    p2 = histogram(theta[burn_in:size(theta)[1], 1], label = "alpha")
    p3 = plot(theta[:,2] .+ No, label = "N")
    p4 = histogram(theta[burn_in:size(theta)[1], 2] .+ No, label = "N")
    plot(p1, p2, p3, p4, layout = (4, 1))
end

N = 10000
a = 0.5
b = (N - 1.0) * a
n = repeat([10, 50], 100)
T = 200
trials = 50
d = truncated(Beta(a, b), upper = 1.0 / maximum(n))
p = rand(d, N)
S = [sample(1:N, pweights(p * i), i, replace = false) for i in n[1:T]]
# S = [sample(1:N, i, replace = false) for i in n[1:T]]
O = Set([i for j in S for i in j])
X = Dict{Any, Any}()
for i in O
    X[i] = sum([i in s for s in S])
end
K = Dict{Int, Int}()
for s in S
    addcounts!(K, s)
end
f = countmap(values(K))
# theta = mcmc(1.0, Float64(length(X)), X, sum(n), 0.0001, 0.0001, 1000; chain = 50000, prop_var_a = 0.5, prop_var_Nu = 0.5)
# res = summarize(theta, length(X), 2500)
# res
# plot_results(theta, length(X), 2500)

@model function demo(X, sum_n)
    Nu ~ LogNormal()
    a ~ Gamma()
    if ((1.0 - exp(loglike(0, sum_n, length(X), a, Nu))) <= 0.0)
        Turing.@addlogprob! -Inf
    else
        Turing.@addlogprob! sum([loglike(X[i], sum_n, length(X), a, Nu) for i in keys(X)]) - length(X) * log(1.0 - exp(loglike(0, sum_n, length(X), a, Nu)))
    end
end

chain = sample(demo(X, sum(n)), NUTS(), 5000)
println(mean(get(chain, [:Nu, :a])[:Nu] .+ length(X)))
histogram(get(chain, [:Nu, :a])[:Nu] .+ length(X))

