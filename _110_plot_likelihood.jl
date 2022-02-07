include("_100_utils.jl")

using .Utils

import Distributions: Beta
import Random: seed!
import Plots

const FIGURES_FOLDER = "./_900_output/figures/"
const N = 100; # True population size
const n = 25; # Sample size at each replication
const T = 2; # Number of replications
# Matrix of true alpha and beta parameters for simulation
const PARAMS = [1.0 4.0; 1.0 1.0; 4.0 1.0];


seed!(555);
Plots.pyplot();
alpha = beta = 0.1:0.5:15; # Alpha and beta values for plotting
for i in 1:size(PARAMS)[1]
    d = Beta(PARAMS[i, 1], PARAMS[i, 2]);
    p = rand(d, N);
    K = simulate_caprecap(N, n, T, p);
    LL(x, y) = loglikelihood(x, y, K, T, "full");
    L(x, y) = likelihood_simple(x, y, K, T, "full");
    flabels = Dict(LL => "Log-likelihood", L => "Likelihood")
    for f in (L, LL)
        title = "$(flabels[f]), true  α = $(PARAMS[i, 1]), β = $(PARAMS[i, 2])";
        Plots.plot(alpha, beta, f,
                st=:contourf, title=title, fontsize=12);
        Plots.plot!(size=(900, 900));
        Plots.xlabel!("α");
        Plots.ylabel!("β");
        Plots.scatter!([PARAMS[i, 1]], [PARAMS[i, 2]],
                    color="lime", label="True α,β");
        Plots.savefig("$(FIGURES_FOLDER)$(flabels[f])_contplt_$i.pdf")
        Plots.plot(alpha, beta, f,
                st=:surface, title=title, fontsize=12,
                camera=(-30, 30), colorbar_title="$(flabels[f])");
        Plots.plot!(size=(1000, 900));
        Plots.xlabel!("α");
        Plots.ylabel!("β");
        Plots.savefig("$(FIGURES_FOLDER)$(flabels[f])_sfplt_$i.pdf")
    end
end
