using .Utils

using Distributions
using SpecialFunctions
using StatsBase
using NLopt

ALPHA = 1;
N = 10000;
T = 2;
AVG_SAMPLE_SIZE = 250;

n = rand(Poisson(AVG_SAMPLE_SIZE), T);
d = Beta(ALPHA, ALPHA * (N - 1));
p = rand(d, N);
S = [sample(1:N, ProbabilityWeights(p), i, replace=false) for i in n];
O = Set([i for j in S for i in j]);
LL(x, grad) = -loglh_truncated(x[1], x[2], S, 100000);
opt = Opt(:LN_SBPLX, 2);
lower = [0.1, 2];
opt.lower_bounds = lower;
opt.ftol_abs = 0.01;
opt.min_objective = LL;
(minf, minx, ret) = NLopt.optimize(opt, [5.0, 10.0]);
