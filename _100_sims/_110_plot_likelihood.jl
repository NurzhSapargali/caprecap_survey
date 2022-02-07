import StatsBase: sample, countmap, Weights, addcounts!
import SpecialFunctions: gamma

N = 100; # True population size
n = 20; # Sample size at each replication
T = 2; # Number of replications


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
        cap = sample(pop, Weights(prob), sample_size, replace=false);
        addcounts!(K, cap);
    end
    return K
end


function likelihood(alpha, beta,
                    freqs::Dict{Int, Int}, reps::Int)
    lh = 1.0;
    for k in values(freqs)
        numerator = ( gamma(alpha + k) 
                    * gamma(beta + reps - k) 
                    * gamma(alpha + beta) );
        denom = ( gamma(alpha) * gamma(beta) * gamma(alpha + beta + reps)
                - gamma(alpha) * gamma(beta + reps) * gamma(alpha + beta) );
        lh *= binomial(reps, k) * numerator / denom;
    end
    return lh
end


function likelihood_simple(alpha, beta,
                           freqs::Dict{Int, Int}, reps::Int)
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
        thirdprod *= ( binomial(reps, k) * multiplier
                     * firstinnerprod * secondinnerprod );
    end
    return (firstprod - secondprod)^(-length(K)) * thirdprod;
end
