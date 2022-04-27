using .Utils

using Distributions
using SpecialFunctions
using StatsBase
using NLopt

import StatsBase: efraimidis_a_wsample_norep!

ALPHAS = [0.1];
N = 1000;
T = [5];
AVG_SAMPLE_SIZE = 50.0;


estis = [];
links = [];
schnab = [];
for alpha in ALPHAS
    d = Beta(alpha, alpha * (N - 1));
    p = rand(d, N);
    for t in T
        n = rand(Poisson(AVG_SAMPLE_SIZE), t);
        for trial in 1:100
            println("***TRIAL NO $trial, $alpha, $t***")
            S = [efraimidis_a_wsample_norep!(1:N, Weights(p), zeros(i)) for i in n];
            O = Set([i for j in S for i in j]);
            LL(x, grad) = -loglh_truncated(x[1], x[2], S, 1000);
            opt = Opt(:LN_SBPLX, 2);
            lower = [0.01, 0];
            opt.lower_bounds = lower;
            opt.min_objective = LL;
            opt.xtol_abs = 0.5;
            (minf, minx, ret) = NLopt.optimize(opt, [5.0, 10.0]);
            push!(estis, [round(minx[2]) + length(O), alpha, t]);
            if t == 2
                push!(links, [round(lincoln(S)), alpha, t]);
            end
            if t > 2
                push!(schnab, [round(schnabel(S)), alpha, t]);
            end
        end
    end
end
