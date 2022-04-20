using .Utils

using Distributions
using SpecialFunctions
using StatsBase
using Plots

ALPHA = 1;
N = 10000;
T = 2;
AVG_SAMPLE_SIZE = 250;

n = rand(Poisson(AVG_SAMPLE_SIZE), T)
d = Beta(ALPHA, ALPHA * (N - 1));
p = rand(d, N);
S = [sample(1:N, ProbabilityWeights(p), i, replace=false) for i in n]
O = Set([i for j in S for i in j]);
LL(alpha, N_u) = -loglh_truncated(alpha, N_u, S, 5000);
xs = [];
lower = [0.1, 2];
upper = [Inf, Inf];
cbk = tr -> begin
            push!(xs, tr[end].metadata["centroid"])
            false
            end
res = optimize(LL, lower, upper, [5.0, 10.0], Fminbox(NelderMead()), Optim.Options(store_trace=true, extended_trace=true, iterations=3, callback=cbk));
