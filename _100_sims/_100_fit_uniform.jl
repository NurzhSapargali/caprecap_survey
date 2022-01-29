import Plots
import Distributions
import StatsBase

N = 1000;
T = 5;
n = 200;


function loglikelihood(K, T, alpha, beta)
    betaT = beta + T;
    firstprod = 1;
    secondprod = 1;
    for j in 1:T
        firstprod *= (1.0 - (j - alpha) / betaT);
        secondprod *= (1.0 - j / betaT);
    end
    firstsum = sum([k * log(alpha + k) for k in values(K)]);
    secondsum = log(betaT) * sum(values(K));
    thirdsum = sum([log(1.0 - j / (alpha + k)) for k in values(K) for j in 1:k]);
    fourthsum = sum([log(1.0 - (k - j) / betaT) for k in values(K) for j in 1:(T-k)]);
    loglh = (length(K) * log(firstprod - secondprod) + firstsum + secondsum
             + thirdsum + fourthsum);
    println(firstprod, ", ", secondprod, ", ",
            firstsum, ", ", secondsum, ", ",
            thirdsum, ", ", fourthsum, ", ",
            loglh);
    return loglh;
end

function main()
    probs = collect(1:N) ./ sum(1:N);
    x = StatsBase.sample(1:N, StatsBase.Weights(probs), n);
    for t in 2:T
        append!(x, StatsBase.sample(1:N, StatsBase.Weights(probs), n));
    end
    k = StatsBase.countmap(x);
    println(loglikelihood(k, T, 1.0, 1.0));
end


main();

