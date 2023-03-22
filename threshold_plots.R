dat_name = "./_900_output/data/diffp/estimates_removal.csv"
dat = read.csv(dat_name)
dat[dat == "Nonsequential frequencies"] = NA
dat[dat == "Noninvertible X'WX"] = NA
dat[dat == Inf] = NA
dat$N_hat = as.numeric(dat$N_hat)
dat$sq_dev = (dat$N_hat - dat$N)^2
dat = drop_na(dat)
sub_dat = dat[(dat$type == "Pseudolikelihood"),]
agg = sub_dat %>%
  group_by(thres) %>%
  summarise(
    avg_log_alpha_hat = mean(log(a_hat)),
    mse = mean(sq_dev),
    hat_mean = mean(N_hat),
    trials = length(N_hat),
    hat_var = var(N_hat)
  )
ggplot(data = agg[2:6,], mapping = aes(x = thres, y = avg_log_alpha_hat)) +
  geom_line() +
  geom_point() +
  geom_hline(aes(yintercept = log(0.5), color = "True log alpha")) +
  geom_hline(aes(yintercept = 0.814, color = "No thresholding")) +
  scale_colour_manual(values = c("blue", "red")) +
  theme_minimal() +
  xlab("Threshold") +
  ylab("") +
  ggtitle("Mean of log alpha estimate, alpha = 0.5, T = 20")
ggsave(filename="mean_alpha_hat.pdf", device = cairo_pdf, width = 297, height = 210, units = "mm")
agg$rmse = sqrt(agg$mse)
agg$bias = agg$hat_mean - 1000
agg$rel_bias = agg$bias / 1000

ggplot(data = agg[2:6,], mapping = aes(x = thres, y = rel_bias)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.0, col = "grey") +
  geom_hline(aes(yintercept = -0.470, color = "No thresholding")) +
  scale_colour_manual(values = c("blue")) +
  theme_minimal() +
  ylab("") +
  xlab("Threshold") +
  ggtitle("Relative bias with thresholding, alpha = 0.5, T = 20")
ggsave(filename="rel_bias_thres.pdf", device = cairo_pdf, width = 297, height = 210, units = "mm")

we_rmse = agg[agg$thres == -999,]$rmse
agg$log_rel_rmse = log(agg$rmse / we_rmse)

ggplot(data = agg[2:6,], mapping = aes(x = thres, y = log_rel_rmse)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0.0, col = "grey") +
  theme_minimal() +
  ylab("") +
  xlab("Threshold") +
  ggtitle("Log of relative RMSE, alpha = 0.5, T = 20")
ggsave(filename="rel_rmse_thres.pdf", device = cairo_pdf, width = 297, height = 210, units = "mm")

dat$type[dat$thres == 3] = "Pseudolikelihood, 3"
dat$type[dat$thres == 4] = "Pseudolikelihood, 4"
dat$type[dat$thres == 5] = "Pseudolikelihood, 5"
dat$type[dat$thres == 6] = "Pseudolikelihood, 6"
dat$type[dat$thres == 7] = "Pseudolikelihood, 7"

agg = dat %>%
  group_by(type) %>%
  summarise(
    mse = mean(sq_dev),
    hat_mean = mean(N_hat),
    trials = length(N_hat),
    hat_var = var(N_hat)
  )
agg$rmse = sqrt(agg$mse)
we_rmse = agg[agg$type == "Pseudolikelihood",]$rmse
agg$log_rel_rmse = log(agg$rmse / we_rmse)
agg$bias = agg$hat_mean - 1000
agg$rel_bias = agg$bias / 1000
agg$type = factor(agg$type, levels = agg$type[order(agg$log_rel_rmse, decreasing = TRUE)])
ggplot(agg[agg$type != "Pseudolikelihood",], mapping = aes(y = type, x = log_rel_rmse)) + 
  geom_bar(stat = "identity") +
  theme_minimal() +
  ggtitle("Log of relative RMSE, alpha = 0.5, T = 20")
ggsave(filename="rel_rmse_all.pdf", device = cairo_pdf, width = 297, height = 210, units = "mm")
agg$type = factor(agg$type, levels = agg$type[order(agg$rel_bias, decreasing = FALSE)])
ggplot(agg, mapping = aes(y = type, x = rel_bias)) + 
  geom_bar(stat = "identity") +
  theme_minimal() +
  ggtitle("Relative bias, alpha = 0.5, T = 20")
ggsave(filename="rel_bias_all.pdf", device = cairo_pdf, width = 297, height = 210, units = "mm")
