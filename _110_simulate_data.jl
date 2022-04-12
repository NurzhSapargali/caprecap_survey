using Distributions
using SpecialFunctions
using Optim
using StatsBase

ALPHA = 2.0;
N = 100;
n = [5; 5; 10; 4; 3; 7; 16; 10; 13; 27; 2; 55; 75; 7; 3; 4; 30; 5; 5; 10; 4; 3; 7; 16; 10; 13; 27; 2; 55; 75; 7; 3; 4; 30]


d = Beta(ALPHA, ALPHA * (N - 1));
p = rand(d, N);
S = [sample(1:N, ProbabilityWeights(p), i, replace=false) for i in n]
