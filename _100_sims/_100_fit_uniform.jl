using Optim, NLSolversBase

import StatsBase
import Random

import Plots
import Distributions

N = 1000;
T = 5;
n = 50;


function loglikelihood(freqs, reps, params)
    betaT = params[2] + reps;
    firstprod = 1;
    secondprod = 1;
    for j in 1:reps
        firstprod *= (1.0 - (j - params[1]) / betaT);
        secondprod *= (1.0 - j / betaT);
    end
    firstsum = log(betaT) * sum(freqs);
    secondsum = 0;
    thirdsum = 0;
    fourthsum = 0;
    for k in freqs
        secondsum += k * log(params[1] + k);
        for j in 1:k
            thirdsum += log(1.0 - j / (params[1] + k));
        end
        for j in 1:(reps - k)
            fourthsum += log(1.0 - (k - j) / betaT);
        end
    end
    loglh = (-length(freqs) * log(firstprod - secondprod) - firstsum + secondsum
             + thirdsum + fourthsum);
    println("alpha = $(params[1]), beta = $(params[2]), loglh = $loglh")
    return -loglh;
end

function main(verbose = true)
    Random.seed!(111);
    x = StatsBase.sample(1:N, n, replace=false);
    for t in 2:T
        append!(x, StatsBase.sample(1:N, n, replace=false));
    end
    K = values(StatsBase.countmap(x));
    println("Optimisation started...");
    f(x) = loglikelihood(K, T, x);
    lower = [1e-8, 1e-8];
    upper = [Inf, Inf];
    opt = optimize(f, lower, upper, x0, NelderMead());
    alpha, beta = Optim.minimizer(opt);
    return
end


main();

