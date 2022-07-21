include("_100_utils.jl")
using .Utils

using Distributions
using SpecialFunctions
using StatsBase
using NLopt
#using Evolutionary

import Random: seed!

N = 3000;
T = [2, 5, 10, 15, 20];
TRIALS = 100;
ALPHAS = [0.1, 3, 10];
AVG_SAMPLE_SIZE = 30;


estis = [];
links = [];
schnab = [];
chaos = [];
chaos_corr = [];
jks = [];
seed!(111);
for alpha in ALPHAS
    d = Beta(alpha, (N - 1.0) * alpha);
    p = rand(d, N);
    n = ones(maximum(T)) * AVG_SAMPLE_SIZE;
    for trial in 1:TRIALS
        println("Generating samples for $n")
        samples = repeat([[]], maximum(T));
        Threads.@threads for i = 1:maximum(T)
           samples[i] = sampford_sample(p, Int.(n[i]));
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
            LL(x, grad) = -loglh_truncated(x[1], x[2], S, O, t, n[1:t], 1000);
            opt = Opt(:LN_SBPLX, 2);
            lower = [0.01, 0];
            opt.lower_bounds = lower;
            opt.min_objective = LL;
            opt.xtol_abs = 0.1;
            (minf, minx, ret) = NLopt.optimize(opt, [5.0, length(O)]);
#             LL(x) = -loglh_truncated(x[1], x[2], S, O, t, n[1:t], 3000);
#             x0 = [1.0, length(O)];
#             lower = [0.01, 0];
#             upper = [Inf, Inf];
#             ga = GA(populationSize=100,selection=uniformranking(3),
#                     mutation=gaussian(),crossover=uniformbin())
#             results = Evolutionary.optimize(LL,
#                                             BoxConstraints(lower, upper),
#                                             x0,
#                                             ga,
#                                             Evolutionary.Options(abstol=1.0, iterations=1000));
#             minx = Evolutionary.minimizer(results);
            push!(estis, [minx[1] ,round(minx[2]) + length(O), t]);
            write_row("_900_output/data/diffp_eqn/estis_$(alpha).csv",
                    [minx[1] ,round(minx[2]) + length(O), t, trial]);
            alpha_trace = [loglh_truncated(i, N - length(O), S, O, t, n[1:t], 1000) for i in 0.1:1:25.1];
            write_row("_900_output/data/diffp_eqn/alpha_trace_$(alpha).csv",
                      vcat(alpha_trace, [minx[1], round(minx[2]), length(O), trial]));
            Nu_trace = [loglh_truncated(alpha, i, S, O, t, n[1:t], 1000) for i in 0:100:5000];
            write_row("_900_output/data/diffp_eqn/Nu_trace_$(alpha).csv",
                      vcat(Nu_trace, [minx[1], round(minx[2]), length(O), trial]));
            push!(chaos, [round(chao(length(O), f)), t, trial]);
            write_row("_900_output/data/diffp_eqn/chaos_$(alpha).csv",
                      [round(chao(length(O), f)), t, trial]);
            push!(chaos_corr, [round(chao_corrected(length(O), t, f)), t, trial]);
            write_row("_900_output/data/diffp_eqn/chaos_corr_$(alpha).csv", 
                      [round(chao_corrected(length(O), t, f)), t, trial]);
            push!(jks, vcat([round(jackknife(length(O), t, f, k)) for k in 1:5], [t, trial]));
            write_row("_900_output/data/diffp_eqn/jks_$(alpha).csv",
                      vcat([round(jackknife(length(O), t, f, k)) for k in 1:5], [t, trial]));
            if t == 2
                push!(links, [round(lincoln(S, n[1:t])), t]);
                write_row("_900_output/data/diffp_eqn/links_$(alpha).csv",
                          [round(lincoln(S, n[1:t])), t, trial]);
            end
            if t > 2
                push!(schnab, [round(schnabel(S, n[1:t])), t]);
                write_row("_900_output/data/diffp_eqn/schnab_$(alpha).csv",
                          [round(schnabel(S, n[1:t])), t, trial]);
            end
        end
    end
end
