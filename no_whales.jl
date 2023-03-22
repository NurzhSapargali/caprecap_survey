include("_110_estimator.jl")
include("_100_utils.jl")
include("_120_benchmarks.jl")

using .Estimator
using .Utils
using .Benchmarks

using StatsBase
using DelimitedFiles

ALPHAS::Float64 = 0.5
DATA_FOLDER::String = "./_200_input/diffp/"
breaks_T::Vector{Int64} = [5, 10, 15, 20]
OUTPUT_FOLDER::String = "./_900_output/data/diffp/"
GRID_SIZE::Int64 = 10000

function get_truth(filename)
    data = readdlm(filename, ' ', String, '\n')[2, 1]
    return parse(Float64, split(data, ",")[1])
end

output_file = OUTPUT_FOLDER * "estimates_removal.csv"
write_row(output_file, ["a_hat", "N_hat", "Nu_hat", "No", "No_cut", "whales", "trial", "T", "alpha", "N", "type", "thres"])
data_folder = DATA_FOLDER * "alpha_0.5/"
data_files = data_files = ["sample_$(i).csv" for i in 1:5000]
N = get_truth(data_folder * "metadata_0.5.csv")
for file in data_files
    trial_no = parse(Int, split(split(file, "_")[2], ".")[1])
    S = read_captures(data_folder * file)
    K = Dict{Int, Int}()
    for s in S
        addcounts!(K, s)
    end
    f = countmap(values(K))
    println("***TRIAL NO $file, 20***")
    O = Set([i for j in S for i in j])
    n = [length(s) for s in S]
    (minf, minx, ret) = fit_model(S, O, n, GRID_SIZE)
    write_row(output_file,
              [minx[1], minx[2] + length(O), minx[2], length(O), 0, 0, trial_no, 20, 0.5, N, "Pseudolikelihood", -999.0])
    for lim in [3, 4, 5, 6, 7]
        whales = [i for i in keys(K) if K[i] >= lim]
        cut_O = symdiff(O, whales)
        cut_n = [length(setdiff(Set(s), Set(whales))) for s in S]
        (minf, minx, ret) = fit_model(S, cut_O, cut_n, GRID_SIZE)
        write_row(output_file,
                  [minx[1], minx[2] + length(cut_O) + length(whales), minx[2], length(O), length(cut_O), length(whales), trial_no, 20, 0.5, N, "Pseudolikelihood", lim])
    end
    benchmarks = Dict{}()
    benchmarks["Schnabel"] = schnabel(S, n)
    benchmarks["Chao"] = chao(length(O), f)
    benchmarks["Zelterman"] = zelterman(length(O), f)
    benchmarks["Conway-Maxwell-Poisson"] = conway_maxwell(length(O), f)
    benchmarks["Huggins"] = huggins(20, K)
    benchmarks["Turing Geometric"] = turing_geometric(length(O), f, 20)
    benchmarks["Turing"] = turing(length(O), f, 20)
    for k in 1:5
        jk = jackknife(length(O), 20, f, k)
        benchmarks["Jackknife k = $(k)"] = jk
    end
    for b in keys(benchmarks)
        write_row(output_file,
                  [-999.0, benchmarks[b], -999.0, length(O), 0, 0, trial_no, 20, 0.5, N, b, -999.0])
    end
end
