module Utils

using StatsBase

export sampford_sample, write_row, freq_of_freq,
       cap_freq, read_captures, lower_pr


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

function pareto_sample(p::Vector{Float64}, n::Int64)
    pr = n * p
    N = length(pr)
    U = rand(N)
    Q = ( U ./ (1.0 .- U) ) ./ ( pr ./ (1.0 .- pr) )
    return sortperm(Q)[1:n]
end

function ar_pareto_sample(p::Vector{Float64}, n::Int64)
    pr = n * p
    N = length(pr)
    d = sum(pr .* (1.0 .- pr))
    sigma2_k = 1.0 ./ (d .+ pr .* (1.0 .- pr))
    ckc0 = (1.0 .- pr) .* sigma2_k.^0.5 .* exp.(sigma2_k .* pr.^2 / 2)
    ckc0 = ckc0 / sum((1.0 .- pr) .* sigma2_k.^0.5 .* exp.(sigma2_k .* pr.^2 / 2))
    ckc0 = (N - n) * ckc0
    Jk = ckc0 ./ (1.0 .- pr)
    A = minimum(Jk)
    S = []
    accept = false
    while (!accept)
        S = pareto_sample(p, n)
        car = sum(1.0 .- pr[S]) / sum(Jk[S] / A .* (1.0 .- pr[S]))
        println(car)
        U = rand()
        if (U <= car)
            accept = true
        end
    end
    return S
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

function lower_pr(N, n, q)
    return exp((n - 1) * ( log(N * (1.0 - q) - 1.0) - log(N) - log(1.0 - q) ) + log(n - 1.0 + N * (1.0 - q)) - log(N) - log(1.0 - q))
end

end
