include("_100_utils.jl")

using .Utils

import Distributions: Beta, pdf
import Random: seed! 

using Optim
using Plots

const FIGURES_FOLDER = "./_900_output/figures/"
const N = 100; # True population size
const n = 10; # Number of replications
const T = 5;
# Matrix of true alpha parameters for simulation
const PARAMS = [0.25, 1.0, 4.0];


seed!(555);
q = (N / n - 1.0);
for alpha in PARAMS
    beta = alpha * q;
    d = Beta(alpha, beta);
    p = rand(d, N);
    p = p ./ sum(p) .* n;
    while ((sum(p .>= 1) > 0) || (sum(p .<= 0) > 0))
        p = rand(d, N);
        p = p ./ sum(p) .* n;
    end
    estis = [];
    for i in 1:100
        K = simulate_caprecap(N, n, T, p);
        println("Data simulated");
        function g!(G, x)
            firstprod = 1.0;
            secondprod = 1.0;
            firstsum = 0.0;
            secondsum = 0.0;
            thirdsum = 0.0;
            fourthsum = 0.0;
            fifthsum = 0.0;
            for j in 1:T
                firstprod *= (1.0 - (j - x[1]) / (x[1] * x[2] + T));
                secondprod *= (1.0 - j / (x[1] * x[2] + T));
                firstsum += (x[2] * j + T) * (1.0 - (j - x[1]) / (x[2] * x[1] + T))^(-1.0);
                secondsum += j * (1.0 - j / (x[1] * x[2] + T))^(-1.0);
                thirdsum += (j - x[1]) * (1.0 - (j - x[1]) / (x[1] * x[2]  + T))^(-1.0);
            end
            for k in values(K)
                firstinnersum = 0.0;
                secondinnersum = 0.0;
                for j in 1:(T-k)
                    firstinnersum += (k + j) / (x[2] * x[1] + T - k - j);
                end
                fourthsum += k - firstinnersum;
                for j in 1:k
                    secondinnersum += j / (x[1] + k - j);
                end
                fifthsum += (x[1] + k)^(-1.0) * (k + secondinnersum);
            end
            G[1] = length(K) / ((x[1] * x[2] + T)^2) * (firstprod * firstsum - x[2] * secondprod * secondsum) / (firstprod - secondprod) + x[2] / (x[2] * x[1] + T) * fourthsum + fifthsum;
            G[2] = length(K) / ((x[1] * x[2] + T)^2) * x[1] * (firstprod * thirdsum - secondprod * secondsum) / (firstprod - secondprod) + x[2] / (x[2] * x[1] + T) * fourthsum;
        end
        L(x) = -loglikelihood(x[1], x[2], K, T);
        lower = [0.01, 0.01]
        upper = [Inf, Inf]
        initial_x = [2.0, 2.0]
        inner_optimizer = GradientDescent()
        println("....Solving")
        results = optimize(L, lower, upper, initial_x, Fminbox(inner_optimizer));
        push!(estis, Optim.minimizer(results));
        println("....Solved! $i");
    end
    
    histogram([i[1] for i in estis],
              xlabel="α_hat", label=false,
              size=(1200, 900), xtickfontsize=22,
              xguidefontsize=24, ytickfontsize=22);
    vline!([i],
           label="True α = $i",
           linewidth=12,
           legendfontsize=22);
    savefig("$(FIGURES_FOLDER)hist_alpha_$i.pdf");
    histogram([i[2] for i in estis],
              xlabel="q", label=false,
              size=(1200, 900), xtickfontsize=22,
              xguidefontsize=24, ytickfontsize=22,
              color="purple");
    vline!([q],
           label="True q = $q",
           linewidth=12,
           legendfontsize=22,
           color="gold");
    savefig("$(FIGURES_FOLDER)hist_q_$i.pdf");
end

