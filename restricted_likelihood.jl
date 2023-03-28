include("_110_estimator.jl")
include("_100_utils.jl")

using .Estimator
using .Utils

using StatsBase
using DelimitedFiles

ALPHA::Float64 = 0.5
DATA_FOLDER::String = "./_200_input/diffp/"
breaks_T::Vector{Int64} = [5, 10, 15, 20]
OUTPUT_FOLDER::String = "./_900_output/data/diffp/"
NU_TRACE_RANGE::Vector{Float64} = collect(0:1.0:500)
GRID_SIZE::Int64 = 10000

function get_truth(filename)
    data = readdlm(filename, ' ', String, '\n')[2, 1]
    return parse(Float64, split(data, ",")[1])
end


output_file = OUTPUT_FOLDER * "estimates_univariate_$(ALPHA).csv"
write_row(output_file, ["a_hat", "N_hat", "Nu_hat", "No", "trial", "T", "alpha", "N", "type"])
data_folder = DATA_FOLDER * "alpha_$(ALPHA)/"
data_files = [file for file in readdir(data_folder) if occursin("sample", file)]
N = get_truth(data_folder * "metadata_$(ALPHA).csv")
grid = LinRange(0.0001, 0.9999, GRID_SIZE)
for file in data_files
    trial_no = parse(Int, split(split(file, "_")[2], ".")[1])
    samples = read_captures(data_folder * file)
    for t in breaks_T
        S = samples[1:t];
        println("***TRIAL NO $file, $t***")
        O = Set([i for j in S for i in j])
        n = [length(s) for s in S]
        X = Dict{Any, Vector{Bool}}()
        for i in O
            X[i] = [i in s for s in S]
        end
        N_o = length(X)
        D = reduce(hcat, [[datalh(p, X[i], n) for p in grid] for i in keys(X)])
        D = hcat(D, [datalh(p, zeros(Bool, length(n)), n) for p in grid])
        (minf, minx, ret) = fit_univariate_model(S, O, n, ALPHA, GRID_SIZE, "alpha")
        write_row(output_file,
                  [-999.0, minx[1] + length(O), minx[1], length(O), trial_no, t, ALPHA, N, "Pseudolikelihood one-parameter"])
        # alpha_trace = [loglh(i, minx[2], S, O, n, 1000) for i in ALPHA_TRACE_RANGE];
        # write_row(OUTPUT_FOLDER * "alpha_trace_$(alpha).csv",
        #           vcat(alpha_trace, [minx[1], minx[2], length(O), t, alpha]));
        println("Evaluating likelihood function at different parameter values...")
        Nu_trace = [loglh(ALPHA, i, N_o, D, grid, false) for i in NU_TRACE_RANGE];
        write_row(OUTPUT_FOLDER * "Nu_trace_$(ALPHA).csv",
                  vcat(Nu_trace, [minx[1], length(O), trial_no, t, ALPHA, N]));
        end
    end
end
