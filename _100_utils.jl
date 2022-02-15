module Utils

import StatsBase: addcounts!
import SpecialFunctions: gamma
import Distributions: Multinomial

export simulate_caprecap, likelihood, likelihood_simple, loglikelihood, sampford

function sampford(prob::Array{Float64, 1}, max_iter=500)
    n = sum(prob);
    n = floor(Int, round(n, digits=0));
    N = length(prob);
    y = prob ./ (1.0 .- prob) ./ sum(prob ./ (1.0 .- prob));
    if sum(y) != 1
        y[length(y)] += 1.0 - sum(y);
    end
    step = 0;
    sb = repeat([2], N);
    a = Multinomial(1, prob ./ sum(prob));
    b = Multinomial(n - 1, y);
    while (sum(sb .<= 1) != N && step <= max_iter)
        sb = rand(a, 1) + rand(b, 1);
        step += 1;
    end
    if (sum(sb .<= 1) == N)
        return sb
    else
        return [];
    end
end


function simulate_caprecap(pop_size::Int, sample_size::Int,
                           reps::Int, prob::Array{Float64, 1})
    if sample_size > pop_size
        error("Sample size is larger than population size")
    end
    if length(prob) != pop_size
        error("Length of probabilities array is not equal to population size")
    end
    if round(sum(prob), digits=0) != sample_size
        prob = prob ./ sum(prob) .* sample_size;
    end
    K = Dict{Int, Int}();
    pop = collect(1:pop_size);
    for t in 1:reps
        mask = [];
        while length(mask) == 0
            mask = sampford(prob);
        end
        s = [i for i in pop if mask[i] == 1];
        addcounts!(K, s);
    end
    return K
end


function likelihood(alpha, q,
                    freqs::Dict{Int, Int}, reps::Int,
                    form::String="proportional")
    if !(form in ("proportional", "full"))
        error("Form of likelihood can be only full or proportional");
    end
    lh = 1.0;
    for k in values(freqs)
        numerator = ( gamma(alpha + k) 
                    * gamma(q * alpha + reps - k) 
                    * gamma(alpha + q * alpha) );
        denom = ( gamma(alpha) * gamma(q * alpha) * gamma(alpha + q * alpha + reps)
                - gamma(alpha) * gamma(q * alpha + reps) * gamma(alpha + q * alpha) );
        lh *= numerator / denom;
        if form == "full"
            lh *= binomial(reps, k);
        end
    end
    return lh
end


function likelihood_simple(alpha, q,
                           freqs::Dict{Int, Int}, reps::Int)
    Ap = 1.0;
    Bp = 1.0;
    thirdprod = 1.0;
    for t in 1:reps
        Ap *= ( 1.0 - (t - alpha) / (q * alpha + reps) );
        Bp *= ( 1.0 - t / (q * alpha + reps) );
    end
    for k in values(freqs)
        multiplier = ( (alpha + k) / (q * alpha + reps) )^k;
        firstinnerprod = 1.0;
        secondinnerprod = 1.0;
        for j in 1:k
            firstinnerprod *= ( 1.0 - j / (alpha + k) );
        end
        for j in 1:(reps - k)
            secondinnerprod *= ( 1.0 - (k + j) / (q * alpha + reps) );
        end
        thirdprod *= ( multiplier * firstinnerprod * secondinnerprod );
    end
    return (Ap - Bp)^(-length(freqs)) * thirdprod;
end


function loglikelihood(alpha, q,
                       freqs::Dict{Int, Int}, reps::Int)
    Ap = 1.0;
    Bp = 1.0;
    secondsum = 0.0;
    for t in 1:reps
        Ap *= ( 1.0 - (t - alpha) / (q * alpha + reps) );
        Bp *= ( 1.0 - t / (q * alpha + reps) );
    end
    for k in values(freqs)
        firstinnersum = 0.0;
        secondinnersum = 0.0;
        for j in 1:k
            firstinnersum += log( alpha + k - j );
        end
        for j in 1:(reps - k)
            secondinnersum += log( q * alpha + reps - k - j );
        end
        secondsum += firstinnersum + secondinnersum
    end
    return ( -reps * length(freqs) * log(q * alpha + reps) - length(freqs) * log( Ap - Bp ) + secondsum);
end

end
