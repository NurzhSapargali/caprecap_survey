library(tidyverse)
library(magrittr)
library(directlabels)
library(ggpubr)
# TO DO: Filter out trials with Infs or completely ditch T values with Infs
OUTPUT_FOLDER = "./_900_output/figures/"

preprocess = function(results){
  res = results[(results$N_hat != Inf)&(results$N_hat >= 0)&(results$N_hat < 100000),]
  res = res[(!is.na(res$N_hat))&(!is.nan(res$N_hat)),]
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

plot_results = function(agg, pop, ylim_bias = c(-3.0, 3.0), ylim_rmse = c(-3.0, 3.0), dense = NA){
  if (is.na(dense)){
    title = paste0("N = ", as.character(pop))
  } else {
    title = paste0("N = ", as.character(pop), ", graph density = ", as.character(dense))
  }
  p1 = ggplot(
    agg[(agg$N == pop)&(agg$type != "Pseudolikelihood"),],
    mapping = aes(x = T, y = lrel_rmse, colour = type, group = type)
  ) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = 0.5) +
    theme_minimal() +
    coord_cartesian(ylim = ylim_rmse, xlim = c(NA, 15.0)) +
    guides(colour = "none") +
    xlab("T") +
    theme_pubr() +
    ylab("LogRelRMSE") + 
    ggtitle(title) +
    geom_dl(
      aes(label = type),
      method = list(
        dl.trans(x = x * 1.025, y = y * 1.0), 
        "last.points",
        "bumpup",
        cex = 0.8
      )
    )
  p2 = ggplot(
    agg[(agg$N == pop)&(agg$type != "Pseudolikelihood"),],
    mapping = aes(x = T, y = lrel_bias, colour = type, group = type)
  ) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, linetype = 2, colour = "black", alpha = 0.5) +
    theme_minimal() +
    coord_cartesian(ylim = ylim_bias, xlim = c(NA, 15.0)) +
    guides(colour = "none") +
    xlab("T") +
    theme_pubr() +
    ylab("LogRelAbsBias") +
    ggtitle("") +
    geom_dl(
      aes(label = type),
      method = list(
        dl.trans(x = x * 1.025, y = y * 1.0), 
        "last.points",
        "bumpup",
        cex = 0.8
      )
    )
  list(p1, p2)
}

for (het in c("0.5", "1.0", "5.0", "10.0", "Inf")){
  dat_name = paste0("./_900_output/data/simulated/estimates_", het, ".csv")
  dat = read.csv(dat_name) %>%
    preprocess()
  agg = aggregate_data(dat)
  agg$lrel_rmse = log(agg$rel_rmse)
  agg$lrel_bias = log(agg$rel_bias)
  agg[agg$type == "Conway-Maxwell-Poisson",]$type = "LCMP"
  ps = c()
  for (N in c(1000, 5000, 10000)){
    ylim_bias = c(-3.0, 3.0)
    ylim_rmse = c(-3.0, 3.0)
    if (het == "0.5"){
      ylim_bias = c(-1.5, 4.5)
      ylim_rmse = c(-1.5, 4.5)
    }
    else if (het == "1.0"){
      if (N == 1000){
        ylim_bias = c(-3.0, 4.0)
        ylim_rmse = c(-3.0, 3.5)
      }
      else if (N == 5000){
        ylim_bias = c(-0.25, 6.0)
        ylim_rmse = c(-2.0, 3.0)
      }
      else{
        ylim_bias = c(-1.0, 4.5)
        ylim_rmse = c(-1.5, 3.0)
      }
    }
    else if (het == "5.0"){
      if (N == 1000){
        ylim_bias = c(-2.0, 6.0)
        ylim_rmse = c(-2.0, 4.5)
      }
      else if (N == 5000){
        ylim_bias = c(-4.0, 7.0)
      }
      else{
        ylim_bias = c(-3.0, 4.5)
      }
    }
    else if (het == "10.0"){
      if (N == 1000){
        ylim_bias = c(-3.5, 10.0)
        ylim_rmse = c(-2.0, 4.0)
      }
      else if (N == 5000){
        ylim_bias = c(-5.0, 6.5)
      }
      else{
        ylim_bias = c(-5.5, 5.5)
      }
    }
    else{
      if (N == 1000){
        ylim_bias = c(-4.0, 3.0)
        ylim_rmse = c(-2.0, 3.0)
      }
      else if (N == 5000){
        ylim_bias = c(-3.5, 2.5)
        ylim_rmse = c(-1.5, 1.5)
      }
      else{
        ylim_bias = c(-5.0, 4.0)
        ylim_rmse = c(-1.5, 2.0)
      }
    }    
    ps = c(ps, plot_results(agg, N, ylim_bias = ylim_bias, ylim_rmse = ylim_rmse))
  }
  ggarrange(plotlist = ps, ncol = 2, nrow = 3)
  ggsave(
    paste0(OUTPUT_FOLDER, "estimates_", het, ".pdf"),
    width = 210,
    height = 297,
    units = "mm"
    )
}

