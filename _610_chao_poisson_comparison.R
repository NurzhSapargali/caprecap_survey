library(tidyverse)
library(ggpubr)

TRIALS <- 1000
SEED <- 777
SAMPLE_SIZES <- seq(100, 600, 15)
POP_SIZE <- 1000
OUTPUT_FOLDER <- "./_900_output/figures/"

srswor_caprecap <- function(n1, n2, N) {
  S <- c(sample(1:N, n1, replace = FALSE), sample(1:N, n2, replace = FALSE))
  K <- table(S)
  f <- c(length(K[K == 1]), length(K[K == 2]))
  return(f)
}

chao <- function(f) {
  return(f[1] + f[2] + f[1]^2 / (2 * f[2]))
}

pseudo_approx <- function(f) {
  return(chao(f) - f[2])
}

dN <- function(f, N) {
  n <- f[1] + 2 * f[2]
  No <- n - f[2]
  return(exp(n / N) * (N - No) - N)
}

dN2 <- function(f, N) {
  n <- f[1] + 2 * f[2]
  No <- n - f[2]
  return(exp(n / N) * (n / N * (No / N - 1) + 1) - 1)
}

pseudo_exact <- function(f, tol = 1e-6) {
  n <- f[1] + 2 * f[2]
  No <- n - f[2]
  converged <- FALSE
  new_N <- No
  while (!converged) {
    old_N <- new_N
    dN <- dN(f, old_N)
    dN2 <- dN2(f, old_N)
    new_N <- old_N - dN / dN2
    converged <- abs(new_N - old_N) < tol
  }
  return(new_N)
}

turing_p <- function(f) {
  return(3 / 2 * f[1] + f[2] + f[1]^2 / (2 * f[2]))
}

simulate_sampling <- function(n1, n2 = NA, N, trials = TRIALS) {
  res <- list()
  if (!anyNA(n2)) {
    for (i in n1) {
      for (j in n2) {
        reps <- t(replicate(trials, srswor_caprecap(i, j, N)))
        reps <- cbind(
          reps,
          apply(reps, 1, chao),
          apply(reps, 1, pseudo_approx),
          apply(reps, 1, pseudo_exact),
          apply(reps, 1, turing_p),
          i,
          j
        )
        res[[paste0(i, j)]] <- reps
      }
    }
  } else {
    for (i in n1) {
      reps <- t(replicate(trials, srswor_caprecap(i, i, N)))
      reps <- cbind(
        reps,
        apply(reps, 1, chao),
        apply(reps, 1, pseudo_approx),
        apply(reps, 1, pseudo_exact),
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
        "Pseudolikelihood (approx.)",
        "Pseudolikelihood (exact)",
        "Turing Poisson",
        "n_1",
        "n_2"
      )
    )
  agg <- res %>%
    group_by(n_1, n_2) %>%
    summarise(
      `Chao` = mean(`Chao`),
      `Pseudolikelihood (approx.)` = mean(`Pseudolikelihood (approx.)`),
      `Pseudolikelihood (exact)` = mean(`Pseudolikelihood (exact)`),
      `Turing Poisson` = mean(`Turing Poisson`),
      .groups = "keep"
    )
  return(agg)
}

set.seed(SEED)
first_fixed <- simulate_sampling(c(100), SAMPLE_SIZES, POP_SIZE)
first_fixed$type <- "'n'[1] == 100"
none_fixed <- simulate_sampling(SAMPLE_SIZES, NA, POP_SIZE)
none_fixed$type <- "'n'[1] == 'n'[2]"

agg <- rbind(first_fixed, none_fixed)
agg$`Chao (approximate difference)` <- (
  agg$`Chao` - agg$`Pseudolikelihood (approx.)`
)
agg$`Turing Poisson (approximate difference)` <- (
  agg$`Turing Poisson` - agg$`Pseudolikelihood (approx.)`
)
agg$`Chao (exact difference)` <- agg$`Chao` - agg$`Pseudolikelihood (exact)`
agg$`Turing Poisson (exact difference)` <- (
  agg$`Turing Poisson` - agg$`Pseudolikelihood (exact)`
)
agg <- pivot_longer(
  agg[, c(-3, -4, -5, -6)],
  cols = c(
    `Chao (approximate difference)`,
    `Turing Poisson (approximate difference)`,
    `Chao (exact difference)`,
    `Turing Poisson (exact difference)`
  ),
  names_to = "method",
  values_to = "difference"
)

agg$rel_cov <- (agg$n_1 + agg$n_2) / POP_SIZE
agg <- agg[agg$rel_cov <= 0.75, ]

ggplot(
  agg,
  mapping = aes(x = (n_1 + n_2) / POP_SIZE, y = difference, colour = type)
) +
  geom_line() +
  geom_point() +
  xlab(expression(("n"[1] + "n"[2]) / N)) +
  ylab("Mean difference") +
  scale_colour_manual(
    values = c("red","blue"),
    labels = expression('n'[1] == 100, 'n'[1] == 'n'[2])
  ) +
  geom_hline(yintercept = 0.0, color = "grey", alpha = 0.9) +
  theme_pubr() +
  facet_wrap(
    vars(method),
    scale = "free_y",
  ) +
  theme(legend.title = element_blank(), legend.position = "right") +
  ggtitle("N = 1000 and 1000 replications")

#ggplot(
#  agg,
#  mapping = aes(x = (n_1 + n_2) / POP_SIZE, y = estimate, colour = method)
#) +
#  geom_line() +
#  geom_point() +
#  geom_hline(yintercept = POP_SIZE, color = "grey", alpha = 0.9) +
#  xlab(expression(("n"[1] + "n"[2]) / N)) +
#  theme_pubr() +
#  facet_wrap(
#    vars(type),
#    scales = "free",
#    labeller = label_parsed
#  ) +
#  theme(legend.title = element_blank())

ggsave(
  paste0(OUTPUT_FOLDER, "appendix/compare_chao_poisson.pdf"),
  width = 297,
  height = 175,
  units = "mm"
)
