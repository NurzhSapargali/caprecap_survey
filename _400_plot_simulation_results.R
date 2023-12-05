library(tidyverse)
library(magrittr)
library(directlabels)

OUTPUT_FOLDER = "./_900_output/figures/"
DRAWS = c(5, 10, 15, 20)

preprocess = function(results){
  res = results[(results$N_hat != Inf)&(results$N_hat >= 0),]
  res$sq_dev = (res$N_hat - res$N)^2
  res$T = as.factor(res$T)
  res
}

aggregate_data = function(pr_results){
  agg = pr_results %>%
    group_by(type, T, N) %>%
    summarise(
      mse = mean(sq_dev),
      hat_mean = mean(N_hat),
      trials = length(N_hat),
      hat_var = var(N_hat)
    )
  agg$rmse = sqrt(agg$mse)
  agg$bias = abs(agg$hat_mean - agg$N)
  nexps = length(agg[agg$type == "Pseudolikelihood",]$rmse)
  repeats = dim(agg)[1] / nexps
  agg$rel_bias = agg$bias / rep(agg[agg$type == "Pseudolikelihood",]$bias, repeats)
  agg$rel_var = agg$hat_var / rep(agg[agg$type == "Pseudolikelihood",]$hat_var, repeats)
  agg$rel_rmse = agg$rmse / rep(agg[agg$type == "Pseudolikelihood",]$rmse, repeats)
  agg
}

plot_alpha_hat = function(pr_results, title_ending, filename, true_alpha = NA){
  header = "Distribution of alpha estimates"
  full_title = paste(header, title_ending, sep = ", ")
  p = ggplot(
    data = pr_results[pr_results$type == "Pseudolikelihood",],
    mapping = aes(y = a_hat, x = T, group = T)
    ) +
    geom_boxplot() +
    stat_summary(fun=mean, colour="darkred", geom="point", shape=18, size=3, show.legend=FALSE) +
    theme_minimal() +
    coord_flip() +
    xlab("T") +
    ylab("Gamma alpha estimate") +
    ggtitle(full_title)
  if (!is.na(true_alpha)){
    p = p + geom_hline(yintercept = true_alpha, colour = "red")
  }
  ggsave(filename, plot = p, device = cairo_pdf, width = 297, height = 210, units = "mm")
}

plot_aggregated_data = function(
    agg,
    y,
#    plot_title,
#    xaxis_breaks,
    ylab,
    xlim = c(NA, NA),
    ylim = c(NA, NA)
    ){
  p = ggplot(agg, mapping = aes_string(x = "T", y = y, colour = "type", group = "type")) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, colour = "black", alpha = 0.25) +
    theme_minimal() +
#    scale_x_continuous(breaks = xaxis_breaks) +
    coord_cartesian(ylim = ylim, xlim = xlim) +
    guides(colour = "none") +
    xlab("T") +
    ylab(ylab)
#    ggtitle(plot_title)
  p
}


dat_name = "./_900_output/data/simulated/estimates_0.5.csv"
dat = read.csv(dat_name) %>%
  preprocess()
agg = aggregate_data(dat)
agg$lrel_rmse = log(agg$rel_rmse)
p1 = plot_aggregated_data(
  agg[(agg$N == 5000)&(agg$type != "Pseudolikelihood"),],
  "lrel_rmse",
  "log(Relative RMSE)",
  ylim = c(NA, 2.5),
  xlim = c(NA, NA)
  ) + 
 # geom_dl(
  #  aes(label = type),
   # method = list(dl.trans(x = x * 1.025, y = y * 1.0),  "last.bumpup",
  #  cex = 0.8
  #  )
#  )
ggsave(
  filename = paste(OUTPUT_FOLDER, "eqp/rel_bias.pdf", sep = ""),
  plot = p1,
  device = cairo_pdf,
  width = 297,
  height = 210,
  units = "mm"
  )
p2 = plot_aggregated_data(
  agg[agg$type != "Pseudolikelihood",],
  "log_rel_rmse",
  "Relative RMSE of alternative estimators, simple random sampling",
  DRAWS,
  "Log of relative RMSE",
  ylim = c(NA, 2.5),
  xlim = c(NA, 24)
) + 
  geom_dl(
    aes(label = type),
    method = list(dl.trans(x = x * 1.025, y = y * 1.0),  "last.bumpup",
                  cex = 0.8
    )
  )
ggsave(
  filename = paste(OUTPUT_FOLDER, "eqp/rel_rmse.pdf", sep = ""),
  plot = p2,
  device = cairo_pdf,
  width = 297,
  height = 210,
  units = "mm"
)
p3 = plot_aggregated_data(
  agg[agg$type != "Pseudolikelihood",],
  "log_rel_var",
  "Relative variance of alternative estimators, simple random sampling",
  DRAWS,
  "Log of relative variance",
  ylim = c(NA, 5),
  xlim = c(NA, 24)
) + 
  geom_dl(
    aes(label = type),
    method = list(dl.trans(x = x * 1.025, y = y * 1.0),  "last.bumpup",
                  cex = 0.8
    )
  )
ggsave(
  filename = paste(OUTPUT_FOLDER, "eqp/rel_var.pdf", sep = ""),
  plot = p3,
  device = cairo_pdf,
  width = 297,
  height = 210,
  units = "mm"
)

