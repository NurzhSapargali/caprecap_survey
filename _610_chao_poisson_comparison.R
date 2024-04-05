library(tidyverse)
library(ggpubr)

TRIALS <- 1000
SEED <- 777
SAMPLE_SIZES <- c(100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)
POP_SIZE <- 1000

srswor_caprecap <- function(n1, n2, N) {
  S <- c(sample(1:N, n1, replace = FALSE), sample(1:N, n2, replace = FALSE))
  K <- table(S)
  f <- c(length(K[K == 1]), length(K[K == 2]))
  return(f)
}

chao <- function(f) {
  return(f[1] + f[2] + f[1]^2 / (2 * f[2]))
}

pseudo <- function(f) {
  return(chao(f) - f[2])
}

turing_p <- function(f) {
  return(3 / 2 * f[1] + f[2] + f[1]^2 / (2 * f[2]))
}

simulate_sampling <- function(n1, n2 = NA, N) {
  res <- list()
  if (!anyNA(n2)) {
    for (i in n1) {
      for (j in n2) {
        reps <- t(replicate(TRIALS, srswor_caprecap(i, j, N)))
        reps <- cbind(
          reps,
          apply(reps, 1, chao),
          apply(reps, 1, pseudo),
          apply(reps, 1, turing_p),
          i,
          j
        )
        res[[paste0(i, j)]] <- reps
      }
    }
  } else {
    for (i in n1) {
      reps <- t(replicate(TRIALS, srswor_caprecap(i, i, N)))
      reps <- cbind(
        reps,
        apply(reps, 1, chao),
        apply(reps, 1, pseudo),
        apply(reps, 1, turing_p),
        i,
        i
      )
      res[[paste0(i, i)]] <- reps
    }
  }
  res <- do.call(rbind, res) %>%
    as.data.frame() %>%
    setNames(
      c(
        "f1",
        "f2",
        "Chao",
        "Pseudolikelihood",
        "Turing Poisson",
        "n_1",
        "n_2"
      )
    )
  agg <- res %>%
    group_by(n_1, n_2) %>%
    summarise(
      `Chao` = mean(`Chao`),
      `Pseudolikelihood` = mean(`Pseudolikelihood`),
      `Turing Poisson` = mean(`Turing Poisson`),
      .groups = "keep"
    )
  return(agg)
}

set.seed(SEED)
first_fixed <- simulate_sampling(c(100), SAMPLE_SIZES, POP_SIZE)
first_fixed$type <- "n1 = 100"
second_fixed <- simulate_sampling(SAMPLE_SIZES, c(100), POP_SIZE)
second_fixed$type <- "n2 = 100"
none_fixed <- simulate_sampling(SAMPLE_SIZES, NA, POP_SIZE)
none_fixed$type <- "n1 = n2"

agg <- rbind(first_fixed, none_fixed, second_fixed) %>%
  pivot_longer(
    cols = c(`Chao`, `Pseudolikelihood`, `Turing Poisson`),
    names_to = "method",
    values_to = "estimate"
  )
agg$type <- factor(agg$type, levels = c("n1 = 100", "n2 = 100", "n1 = n2"))

ggplot(
  agg,
  mapping = aes(x = (n_1 + n_2) / POP_SIZE, y = estimate, colour = method)
) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = POP_SIZE, color = "grey", alpha = 0.7) +
  xlab(expression(("n"[1] + "n"[2]) / N)) +
  theme_pubr() +
  facet_wrap(vars(type), scales = "free_y")

ggsave("check.pdf")
