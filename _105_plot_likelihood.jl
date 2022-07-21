include("_100_utils.jl")
using .Utils

using Plots
using Distributions

import Random: seed!

T = [2, 5, 10, 20];
AVERAGE_SAMPLE_SIZE = 30;
ALPHAS = [0.1, 3.0, 10.0];
N = 3000;
ALPHA_PLOT_RANGE = 0.1:1:25.1;
Nu_PLOT_RANGE = 0:100:5100;

seed!(111);
for alpha in ALPHAS
    d = Beta(alpha, (N - 1.0) * alpha);
    p = rand(d, N);
    n = ones(maximum(T)) * AVG_SAMPLE_SIZE;
    println("Generating samples for $n")
    samples = repeat([[]], maximum(T));
    for i in 1:maximum(T)
        samples[i] = sampford_sample(p, Int.(n[i]));
    end
    for t in T
        S = samples[1:t]
        println("***TRIAL NO $alpha, $trial, $t***")
        O = Set([i for j in S for i in j]);
        L(a, Nu) = loglh_truncated(a, Nu, S, O, t, n[1:t], 1000);
        plot(ALPHA_PLOT_RANGE, Nu_PLOT_RANGE, L, st =:contourf, layout = l);
        xlabel!("Alpha")
        ylabel!("N_u")
