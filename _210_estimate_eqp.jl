include("_110_estimator.jl")
include("_100_utils.jl")
include("_120_benchmarks.jl")

using .Estimator
using .Utils
using .Benchmarks

using StatsBase
using DelimitedFiles

DATA_FOLDER::String = "./_200_input/eqp/"
breaks_T::Vector{Int64} = [5, 10, 15, 20]
OUTPUT_FOLDER::String = "./_900_output/data/eqp/"
GRID_SIZE::Int64 = 10000

function get_truth(filename)
    data = readdlm(filename, ' ', String, '\n')[2, 1]
    return parse(Float64, split(data, ",")[1])
end


data_files = [file for file in readdir(DATA_FOLDER) if occursin("sample", file)]
output_file = OUTPUT_FOLDER * "estimates.csv"
write_row(output_file, ["a_hat", "N_hat", "Nu_hat", "No", "trial", "T", "alpha", "N", "type"])
N = get_truth(DATA_FOLDER * "metadata.csv")
for file in data_files
    samples = read_captures(DATA_FOLDER * file)
    trial_no = parse(Int, split(split(file, "_")[2], ".")[1])
    for t in breaks_T
        S = samples[1:t]
        K = Dict{Int, Int}()
        for s in S
            addcounts!(K, s)
        end
        f = countmap(values(K))
        println("***TRIAL NO $file, $t***")
        O = Set([i for j in S for i in j])
        n = [length(s) for s in S]
        (minf, minx, ret) = fit_model(S, O, n, GRID_SIZE)
        write_row(output_file,
                  [minx[1], minx[2] + length(O), minx[2], length(O), trial_no, t, Inf, N, "Pseudolikelihood"])
#         (minf, minx, ret) = fit_model(S, O, n, MC_DRAWS);
#         write_row(OUTPUT_FOLDER * "estimates.csv",
#                     [minx[1], minx[2], length(O), t]);
#         alpha_trace = [loglh(i, minx[2], S, O, n, 1000) for i in ALPHA_TRACE_RANGE];
#         write_row(OUTPUT_FOLDER * "alpha_trace.csv",
#                     vcat(alpha_trace, [minx[1], minx[2], length(O), t]));
#         Nu_trace = [loglh(minx[1], i, S, O, n, 1000) for i in NU_TRACE_RANGE];
#         write_row(OUTPUT_FOLDER * "Nu_trace.csv",
#                     vcat(Nu_trace, [minx[1], minx[2], length(O), t]));
        benchmarks = Dict{}()
        benchmarks["Schnabel"] = schnabel(S, n)
        benchmarks["Chao"] = chao(length(O), f)
        benchmarks["Zelterman"] = zelterman(length(O), f)
        benchmarks["Conway-Maxwell-Poisson"] = conway_maxwell(length(O), f)
        benchmarks["Huggins"] = huggins(t, K)
        benchmarks["Turing Geometric"] = turing_geometric(length(O), f, t)
        benchmarks["Turing"] = turing(length(O), f, t)
        for k in 1:5
            jk = jackknife(length(O), t, f, k)
            benchmarks["Jackknife k = $(k)"] = jk
        end
        for b in keys(benchmarks)
            write_row(output_file,
                      [-999.0, benchmarks[b], minx[2], length(O), trial_no, t, Inf, N, b])
        end
    end
end
