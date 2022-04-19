using .Utils

using Distributions
using SpecialFunctions
using StatsBase
using Plots

ALPHA = 2;
N = 1000;
n = [50; 50; 100; 40; 30; 30; 100]


d = Beta(ALPHA, ALPHA * (N - 1));
p = rand(d, N);
S = [sample(1:N, ProbabilityWeights(p), i, replace=false) for i in n]
LL(alpha, N_u) = loglh_truncated(alpha, N_u, S, 1000);
O = Set([i for j in S for i in j]);
contourf(1.0:20.0, 2:100:2000, LL)
