include("_110_estimator.jl")
include("_100_utils.jl")
include("_120_benchmarks.jl")

using .Estimator
using .Utils
using .Benchmarks

using StatsBase

import Random: seed!

ALPHAS::Vector{Float64} = [1.0, 5.0, 10.0];
DATA_FOLDER::String = "./_200_input/diffp/";
breaks_T::Vector{Int64} = [5, 10, 15, 20];
OUTPUT_FOLDER::String = "./_900_output/data/diffp/";
const ALPHA_TRACE_RANGE = 0.1:1:50.1;
const NU_TRACE_RANGE = 0.0:200.0:5000.0;
MC_DRAWS::Int64 = 5000;
SEED::Int64 = 111;


seed!(SEED);
for alpha in ALPHAS
    folder = DATA_FOLDER * "alpha_$(alpha)/";
    files = [file for file in readdir(folder) if occursin("sample", file)];
    for file in files
        samples = read_captures(folder * file);
        for t in breaks_T
            S = samples[1:t];
            K = Dict{Int, Int}();
            for s in S
                addcounts!(K, s);
            end
            f = countmap(values(K));
            println("***TRIAL NO $file, $t***");
            O = Set([i for j in S for i in j]);
            n = [length(s) for s in S];
            (minf, minx, ret) = fit_model(S, O, n, MC_DRAWS);
            write_row(OUTPUT_FOLDER * "estimates_$(alpha).csv",
                      [minx[1], minx[2], length(O), t, alpha]);
            alpha_trace = [loglh(i, minx[2], S, O, n, 1000) for i in ALPHA_TRACE_RANGE];
            write_row(OUTPUT_FOLDER * "alpha_trace_$(alpha).csv",
                      vcat(alpha_trace, [minx[1], minx[2], length(O), t, alpha]));
            Nu_trace = [loglh(minx[1], i, S, O, n, 1000) for i in NU_TRACE_RANGE];
            write_row(OUTPUT_FOLDER * "Nu_trace_$(alpha).csv",
                      vcat(Nu_trace, [minx[1], minx[2], length(O), t, alpha]));
            bench_filename = "benchmarks_$(alpha).csv"
            schnab = schnabel(S, n);
            write_row(OUTPUT_FOLDER * bench_filename,
                      [schnab, length(O), t, alpha, "Schnabel"]);
            chao_cor = chao_corrected(length(O), f);
            write_row(OUTPUT_FOLDER * bench_filename,
                      [chao_cor, length(O), t, alpha, "Chao bias-corrected"]);
            zelt = zelterman(length(O), f);
            write_row(OUTPUT_FOLDER * bench_filename,
                      [zelt, length(O), t, alpha, "Zelterman"]);
            cm = conway_maxwell(length(O), f);
            write_row(OUTPUT_FOLDER * bench_filename,
                      [cm, length(O), t, alpha, "Conway-Maxwell-Poisson"]);
            hug = huggins(t, K);
            write_row(OUTPUT_FOLDER * bench_filename,
                      [hug, length(O), t, alpha, "Huggins"]);
            alan = turing(length(O), f, t);
            write_row(OUTPUT_FOLDER * bench_filename,
                      [alan, length(O), t, alpha, "Turing"]);
            for k in 1:5
                jk = jackknife(length(O), t, f, k);
                write_row(OUTPUT_FOLDER * bench_filename,
                          [jk, length(O), t, alpha, "Jackknife k = $(k)"]);
            end
        end
    end
end
