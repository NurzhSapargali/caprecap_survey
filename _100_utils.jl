module Utils

import StatsBase: sample, countmap, ProbabilityWeights, addcounts!
import SpecialFunctions: gamma

export simulate_caprecap, likelihood, likelihood_simple, loglikelihood

function simulate_caprecap(pop_size::Int, sample_size::Int,
                           reps::Int, prob::Array{Float64, 1})
    if sample_size > pop_size
        error("Sample size is larger than population size")
    end
    if length(prob) != pop_size
        error("Length of probabilities array is not equal to population size")
    end
    K = Dict{Int, Int}();
    pop = collect(1:pop_size);
    for t in 1:reps
        cap = sample(pop, ProbabilityWeights(prob), sample_size, replace=false);
        addcounts!(K, cap);
    end
    return K
end


function likelihood(alpha, beta,
                    freqs::Dict{Int, Int}, reps::Int,
                    form::String="proportional")
    if !(form in ("proportional", "full"))
        error("Form of likelihood can be only full or proportional");
    end
    lh = 1.0;
    for k in values(freqs)
        numerator = ( gamma(alpha + k) 
                    * gamma(beta + reps - k) 
                    * gamma(alpha + beta) );
        denom = ( gamma(alpha) * gamma(beta) * gamma(alpha + beta + reps)
                - gamma(alpha) * gamma(beta + reps) * gamma(alpha + beta) );
        lh *= numerator / denom;
        if form == "full"
            lh *= binomial(reps, k);
        end
    end
    return lh
end


function likelihood_simple(alpha, beta,
                           freqs::Dict{Int, Int}, reps::Int,
                           form::String="proportional")
    if !(form in ("proportional", "full"))
        error("Form of likelihood can be only full or proportional");
    end
    firstprod = 1.0;
    secondprod = 1.0;
    thirdprod = 1.0;
    for t in 1:reps
        firstprod *= ( 1.0 - (t - alpha) / (beta + reps) );
        secondprod *= ( 1.0 - t / (beta + reps) );
    end
    for k in values(freqs)
        multiplier = ( (alpha + k) / (beta + reps) )^k;
        firstinnerprod = 1.0;
        secondinnerprod = 1.0;
        for j in 1:k
            firstinnerprod *= ( 1.0 - j / (alpha + k) );
        end
        for j in 1:(reps - k)
            secondinnerprod *= ( 1.0 - (k + j) / (beta + reps) );
        end
        thirdprod *= ( multiplier * firstinnerprod * secondinnerprod );
        if form == "full"
            thirdprod *= binomial(reps, k);
        end
    end
    return (firstprod - secondprod)^(-length(freqs)) * thirdprod;
end


function loglikelihood(alpha, beta,
                       freqs::Dict{Int, Int}, reps::Int,
                       form::String="proportional")
    if !(form in ("proportional", "full"))
        error("Form of likelihood can be only full or proportional");
    end
    firstprod = 1.0;
    secondprod = 1.0;
    firstsum = 0.0;
    secondsum = 0.0;
    bin_num = 0.0;
    for t in 1:reps
        firstprod *= ( 1.0 - (t - alpha) / (beta + reps) );
        secondprod *= ( 1.0 - t / (beta + reps) );
        if form == "full"
            bin_num += log(t);
        end
    end
    for k in values(freqs)
        firstsum += k;
        firstinnersum = 0.0;
        secondinnersum = 0.0;
        for j in 1:k
            firstinnersum += log( 1.0 - j / (alpha + k) );
            if form == "full"
                firstinnersum += -log(j);
            end
        end
        for j in 1:(reps - k)
            secondinnersum += log( 1.0 - (k + j) / (beta + reps) );
            if form == "full"
                secondinnersum += -log(j);
            end
        end
        secondsum += k * log(alpha + k) + firstinnersum + secondinnersum;
    end
    return ( -length(freqs) * log( firstprod - secondprod ) 
            + length(freqs) * bin_num
            - log(beta + reps) * firstsum
            + secondsum);
end

end
