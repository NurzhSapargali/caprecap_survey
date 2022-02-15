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
io = open("optim_data.csv", "w");
for alpha in PARAMS
    beta = alpha * q;
    d = Beta(alpha, beta);
    p = rand(d, N);
    p = p ./ sum(p) .* n;
    while ((sum(p .>= 1) > 0) || (sum(p .<= 0) > 0))
        p = rand(d, N);
        p = p ./ sum(p) .* n;
    end
    for i in 1:100
        write(io, "$alpha,$N,$n,$T");
        K = simulate_caprecap(N, n, T, p);
        println("Data simulated");
        function g!(G, x)
            Ap = 1.0;
            As = 0.0;
            Bp = 1.0;
            Bs = 0.0;
            firstsum = 0.0;
            secondsum = 0.0;
            for j in 1:T
                Ap *= (1.0 - (j - x[1]) / (x[1] * x[2] + T) );
                As += 1.0 / (x[1] + x[1] * x[2] + T - j);
                Bp *= (1.0 - j / (x[1] * x[2] + T) );
                Bs += 1.0 / (x[1] * x[2] + T - j);
            end
            for k in values(K)
                firstinnersum = 0.0;
                secondinnersum = 0.0;
                for j in 1:k
                    firstinnersum += 1.0 / (x[1] + k - j);
                end
                for j in 1:(T - k)
                    secondinnersum += 1.0 / (x[1] * x[2] + T - k - j);
                end
                firstsum += firstinnersum + x[2] * secondinnersum;
                secondsum += secondinnersum;
            end
            G[1] = -(-(1.0 + x[2]) * length(K) * Ap * As / (Ap - Bp) + x[2] * length(K) * Bp * Bs / (Ap - Bp) + firstsum);
            G[2] = -(x[1] * length(K) * (Ap * As - Bp * Bs) / (Ap - Bp) + x[1] * secondsum);
        end
        L(x) = -loglikelihood(x[1], x[2], K, T);
        lower = [0.01, 0.01]
        upper = [Inf, Inf]
        initial_x = [2.0, 2.0]
        inner_optimizer = GradientDescent()
        println("....Solving")
        results = optimize(L, lower, upper, initial_x, Fminbox(inner_optimizer));
        converged = Optim.converged(results);
        write(io, ",$converged");
        row = Optim.minimizer(results);
        for r in row
            write(io, ",$r");
        end
        write(io, "\n");
        println("....Solved! $i");
    end
    
#     histogram([i[1] for i in estis],
#               xlabel="α_hat", label=false,
#               size=(1200, 900), xtickfontsize=22,
#               xguidefontsize=24, ytickfontsize=22);
#     vline!([alpha],
#            label="True α = $alpha",
#            linewidth=12,
#            legendfontsize=22);
#     savefig("$(FIGURES_FOLDER)hist_alpha_$alpha.pdf");
#     histogram([i[2] for i in estis],
#               xlabel="q", label=false,
#               size=(1200, 900), xtickfontsize=22,
#               xguidefontsize=24, ytickfontsize=22,
#               color="purple");
#     vline!([q],
#            label="True q = $q",
#            linewidth=12,
#            legendfontsize=22,
#            color="gold");
#     savefig("$(FIGURES_FOLDER)hist_q_$q.pdf");
end
close(io);
