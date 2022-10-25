module Utils

using StatsBase

export sampford_sample, write_row, freq_of_freq,
       cap_freq, read_captures


function write_row(filename::String, row)
    open(filename, "a") do io
        for (i,j) in enumerate(row)
            if i != length(row)
                print(io, j, ",");
            else
                print(io, j, "\n");
            end
        end
    end
end

function read_captures(filename::String)
    raw = open(f -> read(f, String), filename);
    out = Vector{Int64}[];
    for l in split(raw, '\n')
        if length(l) != 0
            push!(out, parse.(Int, split(l, ",")));
        end
    end
    return out;
end

function sampford_sample(p::Vector{Float64}, n::Int64)
    pr = n * p;
    pr = Dict(enumerate(pr));
    s = Set();
    while length(s) != n
        key = sample(1:length(p), Weights(p));
        val = pop!(pr, key);
        new_pr = collect(values(pr));
        new_pr = new_pr ./ (1.0 .- new_pr);
        new_pr = new_pr ./ sum(new_pr);
        s = Set(sample(collect(keys(pr)), Weights(new_pr), n));
        pr[key] = val;
    end
    #println("Generated")
    return collect(s);
end

function cap_freq(S::Vector{Vector{Int64}})
    K = Dict{Int, Int}();
    for s in S
        addcounts!(K, s);
    end
    return K;
end

function freq_of_freq(K::Dict{Int64, Int64})
    return countmap(values(K));
end

end
