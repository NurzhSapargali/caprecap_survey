library(tidyverse)
library(magrittr)
library(directlabels)


datalist = list()
BA_EDGE_PER_NODE = 1:3

for (i in BA_EDGE_PER_NODE){
  filename = sprintf("./_900_output/data/graphs/ba_%s/nodes_estimates.csv", i)
  dat = read.csv(filename)
  dat[dat == "Nonsequential frequencies"] = NA
  dat[dat == "Noninvertible X'WX"] = NA
  dat[dat == Inf] = NA
  dat$N_hat = as.numeric(dat$N_hat)
  dat$sq_dev = (dat$N_hat - dat$N)^2
  header = "Distribution of alpha estimates, T = 20, 1000 nodes, %s edges"
  bins = 1.0 + 3.322 * log(nrow(dat[dat$a_hat > 0.0,])) 
  ggplot(data = dat[dat$a_hat > 0.0,], mapping = aes(x = log(a_hat))) +
    geom_histogram(size = 0.25, colour = "white", bins = bins) +
    theme_minimal() +
    xlab("Log of alpha estimate") +
    ggtitle(sprintf(header, i * 1000))
  ggsave(
    filename=sprintf("alpha_hat_%s_edges.pdf", i * 1000),
    device = cairo_pdf,
    width = 297,
    height = 210,
    units = "mm"
    )
  cut = drop_na(dat)
  agg = cut %>%
    group_by(type) %>%
    summarise(
      mse = mean(sq_dev),
      hat_mean = mean(N_hat),
      trials = length(N_hat),
      hat_var = var(N_hat)
      )
  agg$rmse = sqrt(agg$mse)
  agg$bias = agg$hat_mean - 1000
  agg$rel_bias = agg$bias / 1000
  agg$edges = i * 1000
  we_rmse = agg[agg$type == "Pseudolikelihood",]$rmse
  agg$log_rel_rmse = log(agg$rmse / we_rmse)
  datalist[[i]] = agg
}

total = do.call(rbind, datalist)
ggplot(
  total[!(total$type %in% c("Pseudolikelihood", "Conway-Maxwell-Poisson")),],
  mapping = aes(x = edges, y = log_rel_rmse, colour = type)
  ) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, colour = "black", alpha = 0.25) +
  theme_minimal() +
  scale_x_continuous(breaks = c(1000, 2000, 3000)) +
  coord_cartesian(ylim = c(NA, NA), xlim = c(775, NA)) +
  geom_dl(
    aes(label = type),
    method = list(dl.trans(x = x * 0.90, y = y * 1.03), "first.bumpup", cex = 0.8)
    ) +
  guides(colour = "none") +
  xlab("Edges") +
  ylab("Log of relative RMSE") +
  ggtitle("Relative RMSE of alternative estimators, T = 20")

ggsave(filename="rel_rmse.pdf", device = cairo_pdf, width = 297, height = 210, units = "mm")

ggplot(
  total[!(total$type %in% c("Conway-Maxwell-Poisson")),],
  mapping = aes(x = edges, y = rel_bias, colour = type)
) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0, colour = "black", alpha = 0.25) +
  theme_minimal() +
  scale_x_continuous(breaks = c(1000, 2000, 3000)) +
  coord_cartesian(ylim = c(NA, NA), xlim = c(775, NA)) +
  geom_dl(
    aes(label = type),
    method = list(dl.trans(x = x * 0.90, y = y * 1.03), "first.bumpup", cex = 0.8)
  ) +
  guides(colour = "none") +
  xlab("Edges") +
  ylab("Relative bias") +
  ggtitle("Relative bias of all estimators, T = 20")

ggsave(filename="rel_bias.pdf", device = cairo_pdf, width = 297, height = 210, units = "mm")
