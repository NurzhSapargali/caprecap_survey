include("_100_utils.jl")
using .Utils

using Distributions
using SpecialFunctions
using StatsBase
using NLopt

import Random: seed!

N = 1000;
T = [2, 4, 10];
AVG_SAMPLE_SIZE = 30;
TRIALS = 100;


estis = [];
links = [];
schnab = [];
chaos = [];
chaos_corr = [];
jks = [];
seed!(777);
n = rand(Poisson(AVG_SAMPLE_SIZE), maximum(T));
for trial in 1:TRIALS
    samples = [sample(1:N, i, replace=false) for i in Int.(n)];
    for t in T
        S = samples[1:t]
        println("***TRIAL NO $trial, $t***")
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
        (minf, minx, ret) = NLopt.optimize(opt, [5.0, 10.0]);
        push!(estis, [minx[1] ,round(minx[2]) + length(O), t]);
        write_row("_900_output/data/eqp_diffn/estis.csv",
                  [minx[1] ,round(minx[2]) + length(O), t, trial]);
        alpha_trace = [loglh_truncated(i, round(minx[2]), S, O, t, n[1:t], 1000) for i in 0.1:1:200.1];
        write_row("_900_output/data/eqp_diffn/alpha_trace.csv",
                  vcat(alpha_trace, [minx[1], round(minx[2]), length(O), trial]));
        Nu_trace = [loglh_truncated(minx[1], i, S, O, t, n[1:t], 1000) for i in 0:100:5000];
        write_row("_900_output/data/eqp_diffn/Nu_trace.csv",
                  vcat(Nu_trace, [minx[1], round(minx[2]), length(O), trial]));
        push!(chaos, [round(chao(length(O), f)), t, trial]);
        write_row("_900_output/data/eqp_diffn/chaos.csv",
                  [round(chao(length(O), f)), t, trial]);
        push!(chaos_corr, [round(chao_corrected(length(O), t, f)), t, trial]);
        write_row("_900_output/data/eqp_diffn/chaos_corr.csv", [round(chao_corrected(length(O), t, f)), t, trial]);
        push!(jks, vcat([round(jackknife(length(O), t, f, k)) for k in 1:5], [t, trial]));
        write_row("_900_output/data/eqp_diffn/jks.csv", vcat([round(jackknife(length(O), t, f, k)) for k in 1:5], [t, trial]));
        if t == 2
            push!(links, [round(lincoln(S, n[1:t])), t]);
            write_row("_900_output/data/eqp_diffn/links.csv", [round(lincoln(S, n[1:t])), t, trial]);
        end
        if t > 2
            push!(schnab, [round(schnabel(S, n[1:t])), t]);
            write_row("_900_output/data/eqp_diffn/schnab.csv", [round(schnabel(S, n[1:t])), t, trial]);
        end
    end
end
