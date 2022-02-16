using DataFrames
using DelimitedFiles
using Plots
using Statistics


data = readdlm("optim_data.csv", ',', '\n');
df = DataFrame(data, ["true_alpha", "N", "n", "T", "converged", "alpha", "q"]);
for alpha in unique(df, "true_alpha")[!,:true_alpha]
    cut = df[(df.true_alpha .== alpha).&(df.converged .== true), :]
    att = nrow(cut);
    histogram(cut[:, "q"], label="q_hat", title="Only $att converged, true alpha = $alpha");
    vline!([9.0], color="red", label="True q");
    vline!([mean(cut[:, "q"])], label="Mean q_hat", color="green")
    vline!([median(cut[:, "q"])], label="Median q_hat", color="purple")
    savefig("hist_q_$(alpha).pdf");
    histogram(log.(cut[:, "alpha"]), label="Log alpha_hat", title="Only $att converged, true alpha = $alpha");
    vline!([log(alpha)], color="red", label="Log true alpha");
    vline!([mean(log.(cut[:, "alpha"]))], label="Mean log alpha_hat", color="green")
    vline!([median(log.(cut[:, "alpha"]))], label="Median log alpha_hat", color="purple")
    savefig("hist_alpha_$(alpha).pdf");
end
