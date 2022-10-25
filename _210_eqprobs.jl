include("_100_utils.jl")

using .Utils

import Distributions: Poisson
import StatsBase: sample

import Random: seed!

N::Int = 1000;
T::Int = 20;
T_min::Int = 5;
TRIALS::Int = 100;
AVG_SAMPLE_SIZE::Int = 30;
SEED::Int = 111;
DATA_FOLDER::String = "./_200_input/eqp/";


seed!(SEED);
metafile = DATA_FOLDER * "metadata.csv";
write_row(metafile, vcat(["N", "T"], ["n_$i" for i in 1:T]));
n = rand(Poisson(AVG_SAMPLE_SIZE), T);
write_row(metafile, vcat([N, T], n));
for trial in 1:TRIALS
    println("Generating samples for $n");
    samples = [sample(1:N, i, replace=false) for i in n];
    O = Set{Int64}();
    while !(length(O) < sum([length(s) for s in samples[1:T_min]]))
        samples = [sample(1:N, i, replace=false) for i in n];
        O = Set([i for j in samples[1:T_min] for i in j]);
        println(length(O) - sum([length(s) for s in samples[1:T_min]]))
    end
    file = DATA_FOLDER * "sample_$(trial).csv";
    for s in samples
        write_row(file, s);
    end
end
