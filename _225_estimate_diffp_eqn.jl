include("_110_estimator.jl")
include("_100_utils.jl")

using .Estimator
using .Utils

import Random: seed!

ALPHAS::Vector{Float64} = [0.1, 0.5, 1.0, 10.0];
DATA_FOLDER::String = "./_200_input/diffp_eqn/";
breaks_T::Vector{Int64} = [5, 10, 15, 20];
OUTPUT_FOLDER::String = "./_900_output/data/diffp_eqn/";
MC_DRAWS::Int64 = 5000;
SEED::Int64 = 111;


seed!(SEED);
for alpha in ALPHAS
    folder = DATA_FOLDER * "alpha_$(alpha)/";
    files = [f for f in readdir(folder) if occursin("sample", f)];
    for f in files
        samples = read_captures(folder * f);
        for t in breaks_T
            S = samples[1:t];
            println("***TRIAL NO $f, $t***");
            O = Set([i for j in S for i in j]);
            n = [length(s) for s in S];
            (minf, minx, ret) = fit_model(S, O, n, MC_DRAWS);
            write_row(OUTPUT_FOLDER * "estimates_$(alpha).csv",
                      [minx[1], minx[2], length(O), t, alpha]);
            alpha_trace = [loglh(i, minx[2], S, O, n, 1000) for i in 0.1:1:50.1];
            write_row(OUTPUT_FOLDER * "alpha_trace_$(alpha).csv",
                      vcat(alpha_trace, [minx[1], minx[2], length(O), t, alpha]));
            Nu_trace = [loglh(minx[1], i, S, O, n, 1000) for i in 0.0:200.0:5000.0];
            write_row(OUTPUT_FOLDER * "Nu_trace_$(alpha).csv",
                      vcat(Nu_trace, [minx[1], minx[2], length(O), t, alpha]));
        end
    end
end
