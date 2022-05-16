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
p = ones(N) / N;
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
        write_row("_900_output/data/estis_eqp_diffn.csv", [minx[1] ,round(minx[2]) + length(O), t]);
        push!(chaos, [round(chao(length(O), f)), t]);
        write_row("_900_output/data/chaos_eqp_diffn.csv", [round(chao(length(O), f)), t]);
        push!(chaos_corr, [round(chao_corrected(length(O), t, f)), t]);
        write_row("_900_output/data/chaos_corr_eqp_diffn.csv", [round(chao_corrected(length(O), t, f)), t]);
        push!(jks, vcat([round(jackknife(length(O), t, f, k)) for k in 1:5], [t]));
        write_row("_900_output/data/jks_eqp_diffn.csv", vcat([round(jackknife(length(O), t, f, k)) for k in 1:5], [t]));
        if t == 2
            push!(links, [round(lincoln(S, n[1:t])), t]);
            write_row("_900_output/data/links_eqp_diffn.csv", [round(lincoln(S, n[1:t])), t]);
        end
        if t > 2
            push!(schnab, [round(schnabel(S, n[1:t])), t]);
            write_row("_900_output/data/schnab_eqp_diffn.csv", [round(schnabel(S, n[1:t])), t]);
        end
    end
end
