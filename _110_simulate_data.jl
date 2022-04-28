using .Utils

using Distributions
using SpecialFunctions
using StatsBase
using NLopt

import StatsBase: efraimidis_a_wsample_norep!
import Random: seed!

ALPHAS = [0.1, 0.5, 1.0, 5.0, 20.0];
N = 1000;
T = [2, 4, 10];
AVG_SAMPLE_SIZE = 75.0;


function write_row(filename, row)
    open(filename, "a") do io
        for i in row
            if i != row[length(row)]
                print(io, i, ",");
            else
                print(io, i, "\n");
            end
        end
    end
end

estis = [];
links = [];
schnab = [];
chaos = [];
chaos_corr = [];
jks = [];
seed!(777);
for alpha in ALPHAS
    d = Beta(alpha, alpha * (N - 1));
    p = rand(d, N);
    for t in T
        n = rand(Poisson(AVG_SAMPLE_SIZE), t);
        for trial in 1:75
            println("***TRIAL NO $trial, $alpha, $t***")
            S = [efraimidis_a_wsample_norep!(1:N, Weights(p), zeros(i)) for i in n];
            O = Set([i for j in S for i in j]);
            K = Dict{Float64, Int}();
            for s in S
                addcounts!(K, s);
            end
            f = countmap(values(K));
            LL(x, grad) = -loglh_truncated(x[1], x[2], S, O, t, n, 1000);
            opt = Opt(:LN_SBPLX, 2);
            lower = [0.01, 0];
            opt.lower_bounds = lower;
            opt.min_objective = LL;
            opt.xtol_abs = 0.1;
            (minf, minx, ret) = NLopt.optimize(opt, [5.0, 10.0]);
            push!(estis, [minx[1] ,round(minx[2]) + length(O), alpha, t]);
            write_row("estis.csv", [minx[1] ,round(minx[2]) + length(O), alpha, t]);
            push!(chaos, [round(chao(length(O), f)), alpha, t]);
            write_row("chaos.csv", [round(chao(length(O), f)), alpha, t]);
            push!(chaos_corr, [round(chao_corrected(length(O), t, f)), alpha, t]);
            write_row("chaos_corr.csv", [round(chao_corrected(length(O), t, f)), alpha, t]);
            push!(jks, vcat([round(jackknife(length(O), t, f, k)) for k in 1:5], [alpha, t]));
            write_row("jks.csv", vcat([round(jackknife(length(O), t, f, k)) for k in 1:5], [alpha, t]));
            if t == 2
                push!(links, [round(lincoln(S, n)), alpha, t]);
                write_row("links.csv", [round(lincoln(S, n)), alpha, t]);
            end
            if t > 2
                push!(schnab, [round(schnabel(n, t, K)), alpha, t]);
                write_row("schnab.csv", [round(schnabel(n, t, K)), alpha, t]);
            end
        end
    end
end
