include("_100_utils.jl")
using .Utils
using RCall
using StatsBase
using Distributions
using DataFrames
using NLopt

import Random: seed!

N = 3000;
T = [5, 10, 15, 20];
TRIALS = 100;
ALPHAS = [20.0];
AVG_SAMPLE_SIZE = 30;


seed!(111);
for alpha in ALPHAS
    d = Beta(alpha, (N - 1.0) * alpha);
    p = rand(d, N);
    n = ones(maximum(T)) * AVG_SAMPLE_SIZE;
    for trial in 1:TRIALS # Do trial = 8, T = 20 case
        println("Generating samples for $n")
        samples = repeat([[]], maximum(T));
        O = Set();
        while !(length(O) < sum([length(s) for s in samples[1:minimum(T)]]))
            Threads.@threads for i = 1:maximum(T)
                samples[i] = sampford_sample(p, Int.(n[i]));
            end
            O = Set([i for j in samples[1:minimum(T)] for i in j]);
        end
        for t in T
            S = samples[1:t]
            println("***TRIAL NO $alpha, $trial, $t***")
            O = Set([i for j in S for i in j]);
            K = Dict{Int, Int}();
            for s in S
                addcounts!(K, s);
            end
            f = countmap(values(K));
            N_o = length(O);
            LL(x, grad) = -loglh(x[1], x[2], S, O, n[1:t], 5000);
            opt = Opt(:LN_SBPLX, 2);
            lower = [0.01, 0];
            opt.lower_bounds = lower;
            opt.min_objective = LL;
            (minf, minx, ret) = NLopt.optimize(opt, [5.0, length(O)]);
            means = minx;
#             noise_points = simulated_annealing(1000, 2.0, N_o,
#                                                5.0, N_o, S,
#                                                O, n[1:t], 500);
#             R"""
#             df <- as.data.frame($noise_points)
#             b <- mgcv::bam(lh ~ s(alpha, N_u, bs="tp"), data=df,
#                         nthreads = 10, discrete = T)
#             f <- function(x){
#                 -predict(b, list(alpha=x[1], N_u=x[2]))
#             }
#             res <- optim(c(5.0, $N_o), f, lower=c(1e-8, 0), method="L-BFGS-B")$par
#             """
#             @rget res
#             means = [mean(noise_points[500:nrow(noise_points),:].alpha),
#                      mean(noise_points[500:nrow(noise_points),:].N_u)];
            write_row("_900_output/data/diffp_eqn/estis_$(alpha).csv",
                      #[means[1], means[2], res[1], res[2], N_o, t, trial]);
                      [means[1], means[2], N_o, t, trial]);
            alpha_trace = [loglh(i, means[2], S, O, n[1:t], 1000) for i in 0.1:1:30.1];
            write_row("_900_output/data/diffp_eqn/alpha_trace_$(alpha).csv",
                    vcat(alpha_trace, [means[1], means[2], N_o, t, trial]));
            Nu_trace = [loglh(means[1], i, S, O, n[1:t], 1000) for i in 0:200:5000];
            write_row("_900_output/data/diffp_eqn/Nu_trace_$(alpha).csv",
                    vcat(Nu_trace, [means[1], means[2], N_o, t, trial]));
            write_row("_900_output/data/diffp_eqn/chaos_$(alpha).csv",
                    [round(chao(N_o, f)), t, trial]);
            write_row("_900_output/data/diffp_eqn/chaos_corr_$(alpha).csv", [round(chao_corrected(N_o, t, f)), N_o, t, trial]);
            write_row("_900_output/data/diffp_eqn/jks_$(alpha).csv", vcat([round(jackknife(N_o, t, f, k)) for k in 1:5], [N_o, t, trial]));
            if t == 2
                write_row("_900_output/data/diffp_eqn/links_$(alpha).csv", [round(lincoln(S, n[1:t])), N_o, t, trial]);
            end
            if t > 2
                write_row("_900_output/data/diffp_eqn/schnab_$(alpha).csv", [round(schnabel(S, n[1:t])), N_o, t, trial]);
            end
        end
    end
end
