library(tidyverse)
library(magrittr)

dat = read.csv("./_900_output/data/graphs/ba_3/nodes_estimates.csv", header = FALSE)
colnames(dat) = c("a_hat", "N_hat", "Nu_hat", "No", "trial", "T", "avg_n", "N", "type")
dat[dat == "Nonsequential frequencies"] = NA
dat[dat == "Noninvertible X'WX"] = NA
dat[dat == Inf] = NA
dat$N_hat = as.numeric(dat$N_hat)
dat$sq_dev = (dat$N_hat - dat$N)^2
cut = drop_na(dat)
agg = cut %>%
  group_by(type) %>%
  summarise(
    mse = mean(sq_dev),
    hat_mean = mean(N_hat),
    trials = length(N_hat),
    hat_var = var(N_hat))
agg$rmse = sqrt(agg$mse)
agg$bias = agg$hat_mean - 1000
agg$rel_bias = abs(agg$bias) / 1000
ggplot(agg[agg$type != "Conway-Maxwell-Poisson",], mapping = aes(y = reorder(type, -rmse), x = rmse)) + geom_bar(stat = "identity")
