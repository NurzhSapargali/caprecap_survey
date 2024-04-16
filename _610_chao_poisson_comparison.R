library(tidyverse)
library(ggpubr)

TRIALS <- 1000
SEED <- 777
SAMPLE_SIZES <- seq(100, 350, 5)
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
first_fixed$type <- "'n'[1] == 100"
none_fixed <- simulate_sampling(SAMPLE_SIZES, NA, POP_SIZE)
none_fixed$type <- "'n'[1] == 'n'[2]"

agg <- rbind(first_fixed, none_fixed)
agg$`Chao (rel.diff)` <- (
  (agg$`Pseudolikelihood` - agg$`Chao`) / POP_SIZE
)
agg$`Turing Poisson (rel.diff)` <- (
  (agg$`Pseudolikelihood` - agg$`Turing Poisson`) / POP_SIZE
)
agg$`Chao (lower)` <- (
  -agg$n_1 / POP_SIZE * agg$n_2 / POP_SIZE
)
agg$`Chao (upper)` <- (
  (agg$n_1 + agg$n_2) / POP_SIZE 
  - 5 / 2 * agg$n_1 / POP_SIZE * agg$n_2 / POP_SIZE
)
agg$`Turing Poisson (lower)` <- (
  -(agg$n_1 + agg$n_2) / (2 * POP_SIZE)
)
agg$`Turing Poisson (upper)` <- (
  (agg$n_1 + agg$n_2) / (POP_SIZE * 2)
  - 3 / 2 * agg$n_1 / POP_SIZE * agg$n_2 / POP_SIZE
)

ggplot(agg, mapping = aes(x = (n_1 + n_2) / POP_SIZE, y = `Chao (rel.diff)`)) +
  geom_line() +
  geom_line(aes(y = `Chao (lower)`), linetype = "dashed") +
  geom_line(aes(y = `Chao (upper)`), linetype = "dashed") +
  xlab(expression(("n"[1] + "n"[2]) / N)) +
  ylab("Mean difference") +
  geom_hline(yintercept = 0.0, color = "grey", alpha = 0.9) +
  theme_pubr() +
  facet_wrap(
    ~ type,
    scale = "free_x",
  ) +
  theme(legend.title = element_blank(), legend.position = "right") +
  ggtitle("N = 1000 and 1000 replications")

p1 <- ggplot(
  agg,
  mapping = aes(x = (n_1 + n_2) / POP_SIZE, y = `Turing Poisson (rel.diff)`)
) +
  geom_line(linewidth = 1.5) +
  geom_line(aes(y = `Turing Poisson (lower)`), linetype = "dashed") +
  geom_line(aes(y = `Turing Poisson (upper)`), linetype = "dashed") +
  geom_ribbon(
    aes(ymin = `Turing Poisson (lower)`, ymax = `Turing Poisson (upper)`),
    fill = "grey",
    alpha = 0.4
  ) +
  xlab(expression(("n"[1] + "n"[2]) / N)) +
  ylab("Mean relative difference") +
  geom_hline(yintercept = 0.0, color = "grey", alpha = 0.9) +
  theme_pubr() +
  facet_wrap(
    ~ type,
    scale = "free",
    labeller = label_parsed
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    strip.text = element_text(size = 14)
  ) +
  ggtitle("Turing Poisson")

p2 <- ggplot(
  agg,
  mapping = aes(x = (n_1 + n_2) / POP_SIZE, y = `Chao (rel.diff)`)
) +
  geom_line(linewidth = 1.5) +
  geom_line(aes(y = `Chao (lower)`), linetype = "dashed") +
  geom_line(aes(y = `Chao (upper)`), linetype = "dashed") +
  geom_ribbon(
    aes(ymin = `Chao (lower)`, ymax = `Chao (upper)`),
    fill = "grey",
    alpha = 0.4
  ) +
  xlab(expression(("n"[1] + "n"[2]) / N)) +
  ylab("Mean relative difference") +
  geom_hline(yintercept = 0.0, color = "grey", alpha = 0.9) +
  theme_pubr() +
  facet_wrap(
    ~ type,
    scale = "free",
    labeller = label_parsed
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    strip.text = element_text(size = 14)
  ) +
  ggtitle("Chao")

p <- ggarrange(p1, p2, ncol = 1, nrow = 2)
annotate_figure(p, top = text_grob("N = 1000 and 1000 replications", size = 16))

ggsave(
  paste0(OUTPUT_FOLDER, "appendix/compare_chao_poisson.pdf"),
  width = 297,
  height = 2 * 175,
  units = "mm"
)
