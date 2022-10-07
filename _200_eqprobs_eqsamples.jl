include("_100_utils.jl")
using .Utils
using RCall
using StatsBase

import Random: seed!

@rlibrary mgcv

N = 3000;
T = [5, 10, 15];
AVG_SAMPLE_SIZE = 30;
TRIALS = 100;


seed!(777);
n = ones(maximum(T)) * AVG_SAMPLE_SIZE;
for trial in 1:TRIALS
    samples = [sample(1:N, i, replace=false) for i in Int.(n)];
    O = Set([i for j in samples[1:minimum(T)] for i in j]);
    while !(length(O) < sum([length(s) for s in samples[1:minimum(T)]]))
        samples = [sample(1:N, i, replace=false) for i in Int.(n)];
        O = Set([i for j in samples[1:minimum(T)] for i in j]);
    end
    for t in T
        S = samples[1:t]
        println("***TRIAL NO $trial, $t***")
        O = Set([i for j in S for i in j]);
        K = Dict{Int, Int}();
        for s in S
            addcounts!(K, s);
        end
        f = countmap(values(K));
        N_o = length(O);
        noise_points = simulated_annealing(1000, 5.0, 2 * N_o,
                                           10.0, N_o, S,
                                           O, n[1:t], 250);
        R"""
        df <- as.data.frame($noise_points)
        b <- mgcv::bam(lh ~ s(alpha, N_u, bs="tp"), data=df,
                       nthreads = 10, discrete = T)
        f <- function(x){
            -predict(b, list(alpha=x[1], N_u=x[2]))
        }
        res <- optim(c(5.0, $N_o), f, lower=c(1e-8, 0), method="L-BFGS-B")$par
        """
        @rget res
        write_row("_900_output/data/eqp_eqn/estis.csv",
                  [mean(noise_points.N_u), res[2], t]);
        alpha_trace = [loglh(i, mean(noise_points.N_u), S, O, n[1:t], 1000) for i in 0.1:1:30.1];
        write_row("_900_output/data/eqp_eqn/alpha_trace.csv",
                  vcat(alpha_trace, [mean(noise_points.alpha), mean(noise_points.N_u), N_o, trial]));
        Nu_trace = [loglh(mean(noise_points.alpha), i, S, O, n[1:t], 1000) for i in 0:200:5000];
        write_row("_900_output/data/eqp_eqn/Nu_trace.csv",
                  vcat(Nu_trace, [mean(noise_points.alpha), mean(noise_points.N_u), N_o, trial]));
        write_row("_900_output/data/eqp_eqn/chaos.csv",
                  [round(chao(N_o, f)), t, trial]);
        write_row("_900_output/data/eqp_eqn/chaos_corr.csv", [round(chao_corrected(N_o, t, f)), t, trial]);
        write_row("_900_output/data/eqp_eqn/jks.csv", vcat([round(jackknife(N_o, t, f, k)) for k in 1:5], [t, trial]));
        if t == 2
            write_row("_900_output/data/eqp_eqn/links.csv", [round(lincoln(S, n[1:t])), t, trial]);
        end
        if t > 2
            write_row("_900_output/data/eqp_eqn/schnab.csv", [round(schnabel(S, n[1:t])), t, trial]);
        end
    end
end
