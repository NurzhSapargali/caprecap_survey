include("_100_utils.jl")
using .Utils

using Plots
using Distributions

import Random: seed!

T = [2, 5, 10, 15, 20];
AVG_SAMPLE_SIZE = 30;
ALPHAS = [0.1, 3.0, 10.0];
N = 3000;

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
        O = Set([i for j in S for i in j]);
        L(a, Nu) = loglh_truncated(a, Nu, S, O, t, n[1:t], 5000);
        plot(0.1:1:(alpha + 2.1),
             0:100:((N - length(O)) + 500),
             L,
             st =:contourf,
             title="Loglikelihood surface for T = $t and alpha = $alpha",
             size=(1024, 640));
        xlabel!("Alpha");
        ylabel!("N_u");
        savefig("contour_$(t)_$(alpha).pdf")
        plot(0.1:1:(alpha + 2.1),
             0:100:((N - length(O)) + 500),
             L,
             st =:surface,
             title="Loglikelihood surface for T = $t and alpha = $alpha",
             size=(1024, 640));
        xlabel!("Alpha");
        ylabel!("N_u");
        savefig("surface_$(t)_$(alpha).pdf")
    end
end
