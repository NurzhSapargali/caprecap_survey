include("_110_estimator.jl")
include("_100_utils.jl")
include("_120_benchmarks.jl")

using .Estimator
using .Utils
using .Benchmarks

using StatsBase

import Random: seed!

DATA_FOLDER::String = "./_200_input/eqp/";
breaks_T::Vector{Int64} = [5, 10, 15, 20];
OUTPUT_FOLDER::String = "./_900_output/data/eqp/";
const ALPHA_TRACE_RANGE = 0.1:1:50.1;
const NU_TRACE_RANGE = 0.0:200.0:5000.0;
MC_DRAWS::Int64 = 5000;
SEED::Int64 = 111;


seed!(SEED);
files = [file for file in readdir(DATA_FOLDER) if occursin("sample", file)];
for file in files
    samples = read_captures(DATA_FOLDER * file);
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
#         (minf, minx, ret) = fit_model(S, O, n, MC_DRAWS);
#         write_row(OUTPUT_FOLDER * "estimates.csv",
#                     [minx[1], minx[2], length(O), t]);
#         alpha_trace = [loglh(i, minx[2], S, O, n, 1000) for i in ALPHA_TRACE_RANGE];
#         write_row(OUTPUT_FOLDER * "alpha_trace.csv",
#                     vcat(alpha_trace, [minx[1], minx[2], length(O), t]));
#         Nu_trace = [loglh(minx[1], i, S, O, n, 1000) for i in NU_TRACE_RANGE];
#         write_row(OUTPUT_FOLDER * "Nu_trace.csv",
#                     vcat(Nu_trace, [minx[1], minx[2], length(O), t]));
        bench_filename = "benchmarks.csv"
        schnab = schnabel(S, n);
        write_row(OUTPUT_FOLDER * bench_filename,
                  [schnab, length(O), t, "Schnabel"]);
        chao_est = chao(length(O), f);
        write_row(OUTPUT_FOLDER * bench_filename,
                  [chao_est, length(O), t, "Chao"]);
        zelt = zelterman(length(O), f);
        write_row(OUTPUT_FOLDER * bench_filename,
                    [zelt, length(O), t, "Zelterman"]);
        cm = conway_maxwell(length(O), f);
        write_row(OUTPUT_FOLDER * bench_filename,
                    [cm, length(O), t, "Conway-Maxwell-Poisson"]);
        hug = huggins(t, K);
        write_row(OUTPUT_FOLDER * bench_filename,
                    [hug, length(O), t, "Huggins"]);
        alan_geo = turing_geometric(length(O), f, t);
        write_row(OUTPUT_FOLDER * bench_filename,
                    [alan_geo, length(O), t, "Turing Geometric"]);
        alan = turing(length(O), f, t);
        write_row(OUTPUT_FOLDER * bench_filename,
                  [alan, length(O), t, "Turing"]);
        for k in 1:5
            jk = jackknife(length(O), t, f, k);
            write_row(OUTPUT_FOLDER * bench_filename,
                      [jk, length(O), t, "Jackknife k = $(k)"]);
        end
    end
end
