library(tidyverse)
library(magrittr)
library(directlabels)

OUTPUT_FOLDER = "./_900_output/figures/"
DRAWS = c(5, 10, 15, 20)

preprocess = function(results, drop_alpha = FALSE){
  results[results == "Nonsequential frequencies"] = NA
  results[results == "Noninvertible X'WX"] = NA
  results[results == Inf] = NA
  results$N_hat = as.numeric(results$N_hat)
  results$sq_dev = (results$N_hat - results$N)^2
  if (drop_alpha){
    results = select(results, -alpha)
  }
  results
}

aggregate_data = function(pr_results, true_N){
  agg = pr_results %>%
    group_by(type, T) %>%
    summarise(
      mse = mean(sq_dev),
      hat_mean = mean(N_hat),
      trials = length(N_hat),
      hat_var = var(N_hat)
    )
  agg$rmse = sqrt(agg$mse)
  agg$bias = agg$hat_mean - true_N
  agg$rel_bias = agg$bias / true_N
  we_rmse = agg[agg$type == "Pseudolikelihood",]$rmse
  we_var = agg[agg$type == "Pseudolikelihood",]$hat_var
  agg$log_rel_rmse = log(agg$rmse / we_rmse)
  agg$log_rel_var = log(agg$hat_var / we_var)
  agg
}

plot_alpha_hat = function(pr_results, title_ending, filename, true_alpha = NA){
  header = "Distribution of alpha estimates"
  full_title = paste(header, title_ending, sep = ", ")
  p = ggplot(
    data = pr_results[pr_results$type == "Pseudolikelihood",],
    mapping = aes(x = log(a_hat), y = T, group = T)
    ) +
    geom_boxplot() +
    theme_minimal() +
    xlab("Log of alpha estimate") +
    ggtitle(full_title)
  if (!is.na(true_alpha)){
    p = p + geom_vline(xintercept = log(true_alpha), colour = "red")
  }
  ggsave(filename, plot = p, device = cairo_pdf, width = 297, height = 210, units = "mm")
}

plot_aggregated_data = function(
    agg,
    y,
    plot_title,
    xaxis_breaks,
    ylab,
    xlim = c(NA, NA),
    ylim = c(NA, NA)
    ){
  p = ggplot(agg, mapping = aes_string(x = "T", y = y, colour = "type")) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, colour = "black", alpha = 0.25) +
    theme_minimal() +
    scale_x_continuous(breaks = xaxis_breaks) +
    coord_cartesian(ylim = ylim, xlim = xlim) +
    guides(colour = "none") +
    xlab("T") +
    ylab(ylab) +
    ggtitle(plot_title)
  p
}


dat_name = "./_900_output/data/eqp/estimates.csv"
dat = read.csv(dat_name) %>%
  preprocess(drop_alpha = TRUE) %>%
  drop_na()
dat %>% plot_alpha_hat(
  "simple random sampling",
  paste(OUTPUT_FOLDER, "eqp/alpha_estimates.pdf", sep = "")
  )
agg = aggregate_data(dat, 1000)
p1 = plot_aggregated_data(
  agg,
  "rel_bias",
  "Relative bias of estimators, simple random sampling",
  DRAWS,
  "Relative bias",
  ylim = c(NA, 0.5),
  xlim = c(NA, 24)
  ) + 
  geom_dl(
    aes(label = type),
    method = list(dl.trans(x = x * 1.025, y = y * 1.0),  "last.bumpup",
    cex = 0.8
    )
  )
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


dat_name = "./_900_output/data/diffp/estimates_0.5.csv"
dat = read.csv(dat_name) %>%
  preprocess() %>%
  drop_na()
dat %>% plot_alpha_hat(
  "alpha = 0.5",
  paste(OUTPUT_FOLDER, "diffp/alpha_estimates_0.5.pdf", sep = ""),
  true_alpha = 0.5
)
agg = aggregate_data(dat, 1000)
p1 = plot_aggregated_data(
  agg,
  "rel_bias",
  "Relative bias of estimators, alpha = 0.5",
  DRAWS,
  "Relative bias",
  ylim = c(NA, 0.20),
  xlim = c(NA, 24)
) + 
  geom_dl(
    aes(label = type),
    method = list(dl.trans(x = x * 1.025, y = y * 1.0),  "last.bumpup", cex = 0.8)
  )
ggsave(
  filename = paste(OUTPUT_FOLDER, "diffp/rel_bias_0.5.pdf", sep = ""),
  plot = p1,
  device = cairo_pdf,
  width = 297,
  height = 210,
  units = "mm"
)
p2 = plot_aggregated_data(
  agg[agg$type != "Pseudolikelihood",],
  "log_rel_rmse",
  "Relative RMSE of alternative estimators, alpha = 0.5",
  DRAWS,
  "Log of relative RMSE",
  ylim = c(NA, 0.5),
  xlim = c(NA, 24)
) + 
  geom_dl(
    aes(label = type),
    method = list(dl.trans(x = x * 1.025, y = y * 1.0),  "last.bumpup",
                  cex = 0.8
    )
  )
ggsave(
  filename = paste(OUTPUT_FOLDER, "diffp/rel_rmse_0.5.pdf", sep = ""),
  plot = p2,
  device = cairo_pdf,
  width = 297,
  height = 210,
  units = "mm"
)
p3 = plot_aggregated_data(
  agg[agg$type != "Pseudolikelihood",],
  "log_rel_var",
  "Relative variance of alternative estimators, alpha = 0.5",
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
  filename = paste(OUTPUT_FOLDER, "diffp/rel_var_0.5.pdf", sep = ""),
  plot = p3,
  device = cairo_pdf,
  width = 297,
  height = 210,
  units = "mm"
)

dat_name = "./_900_output/data/diffp/estimates_1.0.csv"
dat = read.csv(dat_name) %>%
  preprocess() %>%
  drop_na()
dat %>% plot_alpha_hat(
  "alpha = 1",
  paste(OUTPUT_FOLDER, "diffp/alpha_estimates_1.0.pdf", sep = ""),
  true_alpha = 1
)
agg = aggregate_data(dat, 1000)
p1 = plot_aggregated_data(
  agg,
  "rel_bias",
  "Relative bias of estimators, alpha = 1",
  DRAWS,
  "Relative bias",
  ylim = c(NA, 0.20),
  xlim = c(NA, 24)
) + 
  geom_dl(
    aes(label = type),
    method = list(dl.trans(x = x * 1.025, y = y * 1.0),  "last.bumpup", cex = 0.8)
  )
ggsave(
  filename = paste(OUTPUT_FOLDER, "diffp/rel_bias_1.0.pdf", sep = ""),
  plot = p1,
  device = cairo_pdf,
  width = 297,
  height = 210,
  units = "mm"
)
p2 = plot_aggregated_data(
  agg[agg$type != "Pseudolikelihood",],
  "log_rel_rmse",
  "Relative RMSE of alternative estimators, alpha = 1",
  DRAWS,
  "Log of relative RMSE",
  ylim = c(NA, 0.5),
  xlim = c(NA, 24)
) + 
  geom_dl(
    aes(label = type),
    method = list(dl.trans(x = x * 1.025, y = y * 1.0),  "last.bumpup",
                  cex = 0.8
    )
  )
ggsave(
  filename = paste(OUTPUT_FOLDER, "diffp/rel_rmse_1.0.pdf", sep = ""),
  plot = p2,
  device = cairo_pdf,
  width = 297,
  height = 210,
  units = "mm"
)
p3 = plot_aggregated_data(
  agg[agg$type != "Pseudolikelihood",],
  "log_rel_var",
  "Relative variance of alternative estimators, alpha = 1",
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
  filename = paste(OUTPUT_FOLDER, "diffp/rel_var_1.0.pdf", sep = ""),
  plot = p3,
  device = cairo_pdf,
  width = 297,
  height = 210,
  units = "mm"
)

dat_name = "./_900_output/data/diffp/estimates_5.0.csv"
dat = read.csv(dat_name) %>%
  preprocess() %>%
  drop_na()
dat %>% plot_alpha_hat(
  "alpha = 5",
  paste(OUTPUT_FOLDER, "diffp/alpha_estimates_5.0.pdf", sep = ""),
  true_alpha = 5
)
agg = aggregate_data(dat, 1000)
p1 = plot_aggregated_data(
  agg,
  "rel_bias",
  "Relative bias of estimators, alpha = 5",
  DRAWS,
  "Relative bias",
  ylim = c(NA, 0.20),
  xlim = c(NA, 24)
) + 
  geom_dl(
    aes(label = type),
    method = list(dl.trans(x = x * 1.025, y = y * 1.0),  "last.bumpup", cex = 0.8)
  )
ggsave(
  filename = paste(OUTPUT_FOLDER, "diffp/rel_bias_5.0.pdf", sep = ""),
  plot = p1,
  device = cairo_pdf,
  width = 297,
  height = 210,
  units = "mm"
)
p2 = plot_aggregated_data(
  agg[agg$type != "Pseudolikelihood",],
  "log_rel_rmse",
  "Relative RMSE of alternative estimators, alpha = 5",
  DRAWS,
  "Log of relative RMSE",
  ylim = c(NA, 2),
  xlim = c(NA, 24)
) + 
  geom_dl(
    aes(label = type),
    method = list(dl.trans(x = x * 1.025, y = y * 1.0),  "last.bumpup",
                  cex = 0.8
    )
  )
ggsave(
  filename = paste(OUTPUT_FOLDER, "diffp/rel_rmse_5.0.pdf", sep = ""),
  plot = p2,
  device = cairo_pdf,
  width = 297,
  height = 210,
  units = "mm"
)
p3 = plot_aggregated_data(
  agg[agg$type != "Pseudolikelihood",],
  "log_rel_var",
  "Relative variance of alternative estimators, alpha = 5",
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
  filename = paste(OUTPUT_FOLDER, "diffp/rel_var_5.0.pdf", sep = ""),
  plot = p3,
  device = cairo_pdf,
  width = 297,
  height = 210,
  units = "mm"
)

dat_name = "./_900_output/data/diffp/estimates_10.0.csv"
dat = read.csv(dat_name) %>%
  preprocess() %>%
  drop_na()
dat %>% plot_alpha_hat(
  "alpha = 10",
  paste(OUTPUT_FOLDER, "diffp/alpha_estimates_10.0.pdf", sep = ""),
  true_alpha = 10
)
agg = aggregate_data(dat, 1000)
p1 = plot_aggregated_data(
  agg,
  "rel_bias",
  "Relative bias of estimators, alpha = 10",
  DRAWS,
  "Relative bias",
  ylim = c(NA, 0.20),
  xlim = c(NA, 24)
) + 
  geom_dl(
    aes(label = type),
    method = list(dl.trans(x = x * 1.025, y = y * 1.0),  "last.bumpup", cex = 0.8)
  )
ggsave(
  filename = paste(OUTPUT_FOLDER, "diffp/rel_bias_10.0.pdf", sep = ""),
  plot = p1,
  device = cairo_pdf,
  width = 297,
  height = 210,
  units = "mm"
)
p2 = plot_aggregated_data(
  agg[agg$type != "Pseudolikelihood",],
  "log_rel_rmse",
  "Relative RMSE of alternative estimators, alpha = 10",
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
  filename = paste(OUTPUT_FOLDER, "diffp/rel_rmse_10.0.pdf", sep = ""),
  plot = p2,
  device = cairo_pdf,
  width = 297,
  height = 210,
  units = "mm"
)
p3 = plot_aggregated_data(
  agg[agg$type != "Pseudolikelihood",],
  "log_rel_var",
  "Relative variance of alternative estimators, alpha = 10",
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
  filename = paste(OUTPUT_FOLDER, "diffp/rel_var_10.0.pdf", sep = ""),
  plot = p3,
  device = cairo_pdf,
  width = 297,
  height = 210,
  units = "mm"
)
