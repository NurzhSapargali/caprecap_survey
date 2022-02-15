include("_100_utils.jl")

using .Utils

import Distributions: Beta
import Random: seed!
import Plots

const FIGURES_FOLDER = "./_900_output/figures/"
const N = 100; # True population size
const n = 10; # Sample size at each replication
const T = 5; # Number of replications
# Matrix of true alpha parameters for simulation
const PARAMS = [0.25, 1.0, 4.0];


seed!(555);
Plots.pyplot();
alpha = coef = 0.1:0.5:15.1; # Alpha and q values for plotting
q = (N / n - 1.0);
for i in PARAMS
    d = Beta(i, q * i);
    p = rand(d, N);
    p = p ./ sum(p) .* n;
    while ((sum(p .>= 1) > 0) || (sum(p .<= 0) > 0))
        p = rand(d, N);
        p = p ./ sum(p) .* n;
    end
    K = simulate_caprecap(N, n, T, p);
    LL(x, y) = loglikelihood(x, y, K, T);
    L(x, y) = likelihood_simple(x, y, K, T);
    flabels = Dict(LL => "Log-likelihood", L => "Likelihood")
    for f in (L, LL)
        title = "α = $i, q = $(N / n - 1.0)";
        Plots.plot(alpha, coef, f,
                   st=:contourf, title=title, titlefontsize=24,
                   xlabel="α", xtickfontsize=22, xguidefontsize=24,
                   ylabel="q", ytickfontsize=22, yguidefontsize=24,
                   legendfontsize=18, colorbar_tickfontsize=22,size=(1000, 900));
        Plots.scatter!([i], [q],
                        color="lime", label="True α,q",
                        markersize=12);
        Plots.savefig("$(FIGURES_FOLDER)$(flabels[f])_contplt_$i.pdf")
        Plots.plot(alpha, coef, f,
                   st=:surface, title=title, camera=(-30, 30),
                   colorbar_title="$(flabels[f])", titlefontsize=21,
                   xlabel="α", xtickfontsize=18, xguidefontsize=22,
                   ylabel="q", ytickfontsize=18, yguidefontsize=22,
                   ztickfontsize=12, colorbar_titlefontsize=14, colorbar_tickfontsize=12,
                   size=(1000, 900));
        Plots.savefig("$(FIGURES_FOLDER)$(flabels[f])_sfplt_$i.pdf")
    end
end
