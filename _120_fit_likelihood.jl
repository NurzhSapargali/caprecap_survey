include("_100_utils.jl")

using .Utils

import Distributions: Beta, pdf
import Random: seed! 

using Optim
using Plots

const FIGURES_FOLDER = "./_900_output/figures/"
const N = 100; # True population size
const n = 25; # Number of replications
const T = 2;
# Matrix of true alpha and beta parameters for simulation
const PARAMS = [1.0 4.0; 1.0 1.0; 4.0 1.0];


seed!(555);
for row in 1:size(PARAMS)[1]
    alpha = PARAMS[row, 1];
    beta = PARAMS[row, 2]
    d = Beta(alpha, beta);
    p = rand(d, N);
    estis = [];
    for i in 1:100
        K = simulate_caprecap(N, n, T, p);
        function g!(G, x)
            firstprod = 1.0;
            secondprod = 1.0;
            firstsum = 0.0;
            secondsum = 0.0;
            thirdsum = 0.0;
            fourthsum = 0.0;
            fifthsum = 0.0;
            for j in 1:T
                firstprod *= ( 1.0 - (j - x[1]) / (x[2] + T) );
                secondprod *= ( 1.0 - j / (x[2] + T) );
                firstsum += ( 1.0 - (j - x[1]) / (x[2] + T) )^(-1);
                secondsum += (j - x[1]) * ( 1.0 - (j - x[1]) / (x[2] + T) )^(-1);
                thirdsum += j * ( 1.0 - (j - x[1]) / (x[2] + T) )^(-1);
            end
            for k in values(K)
                firstinnersum = 0.0;
                secondinnersum = 0.0;
                for j in 1:k
                    firstinnersum += ( j / (x[1] + k - j) );
                end
                for j in 1:(T - k)
                    secondinnersum += ( (k - j) / (x[2] + T - k + j) );
                end
                fourthsum += (x[1] + k)^(-1) * (k - firstinnersum);
                fifthsum += k - secondinnersum;
            end
            G[1] = -( -length(K) / (x[2] + T)
                    * firstprod * firstsum
                    / (firstprod - secondprod) + fourthsum );
            G[2] = -( -length(K) / ((x[2] + T)^2)
                    * ((firstprod * secondsum - secondprod * thirdsum)
                    / (firstprod - secondprod)) - (x[2] + T)^(-1) * fifthsum ); 
        end
        L(x) = -loglikelihood(x[1], x[2], K, T);
        lower = [0.01, 0.01]
        upper = [Inf, Inf]
        initial_x = [2.0, 2.0]
        inner_optimizer = GradientDescent()
        results = optimize(L, g!, lower, upper, initial_x, Fminbox(inner_optimizer));
        push!(estis, Optim.minimizer(results));
        println("$row, $i");
    end
    
    histogram([i[1] for i in estis],
              xlabel="α_hat", label=false,
              size=(1200, 900), xtickfontsize=22,
              xguidefontsize=24, ytickfontsize=22);
    vline!([alpha],
           label="True α = $alpha",
           linewidth=12,
           legendfontsize=22);
    savefig("$(FIGURES_FOLDER)hist_alpha_$row.pdf");
    histogram([i[2] for i in estis],
              xlabel="β_hat", label=false,
              size=(1200, 900), xtickfontsize=22,
              xguidefontsize=24, ytickfontsize=22,
              color="purple");
    vline!([beta],
           label="True β = $beta",
           linewidth=12,
           legendfontsize=22,
           color="gold");
    savefig("$(FIGURES_FOLDER)hist_beta_$row.pdf");
end

