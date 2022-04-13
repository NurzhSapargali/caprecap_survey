using Distributions
using SpecialFunctions
using Optim
using StatsBase

ALPHA = 10;
N = 1000;
n = [50; 50; 100; 40; 30; 30; 100; 50; 20; 30; 5; 10; 25; 35; 25; 10; 10; 25; 10; 5; 50; 50; 100; 40; 30; 30; 100; 50; 20; 30; 5; 10; 25; 35; 25; 10; 10; 25; 10; 5]


d = Beta(ALPHA, ALPHA * (N - 1));
p = rand(d, N);
S = [sample(1:N, ProbabilityWeights(p), i, replace=false) for i in n]
LL(alpha, N_u) = loglikelihood(alpha, N_u, S);
O = Set([i for j in S for i in j]);
contourf(1.0:20.0, 2:100:2000, LL)
