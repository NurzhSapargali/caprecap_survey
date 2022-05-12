
include("_100_utils.jl")
using .Utils

using Distributions
using SpecialFunctions
using StatsBase
using NLopt

import Random: seed!

ALPHAS = [0.1];
N = 1000;
T = [2, 4, 10];
AVG_SAMPLE_SIZE = 30;


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
    p = rand(Beta(alpha, alpha * (N - 1)), N);
    nmax = Int(floor(1.0 / maximum(p)));
    n = ones(maximum(T)) .* nmax;
    for trial in 1:100
        samples = [pareto_sampling(p, i) for i in Int.(n)];
        for t in T
            S = samples[1:t]
            println("***TRIAL NO $trial, $alpha, $t***")
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
            push!(estis, [minx[1] ,round(minx[2]) + length(O), alpha, t]);
            write_row("_900_output/data/estis_diffp_eqn_$(alpha)_$(nmax).csv", [minx[1] ,round(minx[2]) + length(O), alpha, t]);
            push!(chaos, [round(chao(length(O), f)), alpha, t]);
            write_row("_900_output/data/chaos_diffp_eqn_$(alpha)_$(nmax).csv", [round(chao(length(O), f)), alpha, t]);
            push!(chaos_corr, [round(chao_corrected(length(O), t, f)), alpha, t]);
            write_row("_900_output/data/chaos_corr_diffp_eqn_$(alpha)_$(nmax).csv", [round(chao_corrected(length(O), t, f)), alpha, t]);
            push!(jks, vcat([round(jackknife(length(O), t, f, k)) for k in 1:5], [alpha, t]));
            write_row("_900_output/data/jks_diffp_eqn_$(alpha)_$(nmax).csv", vcat([round(jackknife(length(O), t, f, k)) for k in 1:5], [alpha, t]));
            if t == 2
                push!(links, [round(lincoln(S, n[1:t])), alpha, t]);
                write_row("_900_output/data/links_diffp_eqn_$(alpha)_$(nmax).csv", [round(lincoln(S, n[1:t])), alpha, t]);
            end
            if t > 2
                push!(schnab, [round(schnabel(S, n[1:t])), alpha, t]);
                write_row("_900_output/data/schnab_diffp_eqn_$(alpha)_$(nmax).csv", [round(schnabel(S, n[1:t])), alpha, t]);
            end
        end
    end
end
