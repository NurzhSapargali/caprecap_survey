library(ggplot2)
library(magrittr)
library(dplyr)

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

turing_b <- function(f) {
  return(f[1] + f[2] + f[1]^2 / (4 * f[2]))
}

turing_p <- function(f) {
  return(3 / 2 * f[1] + f[2] + f[1]^2 / (2 * f[2]))
}

set.seed(SEED)
res <- list()
for (n in SAMPLE_SIZES) {
  reps <- t(replicate(TRIALS, srswor_caprecap(n, n, POP_SIZE)))
  reps <- cbind(
    reps,
    apply(reps, 1, chao),
    apply(reps, 1, pseudo),
    apply(reps, 1, turing_p),
    apply(reps, 1, turing_b),
    n
  )
  res[[as.character(n)]] <- reps
}
res <- do.call(rbind, res) %>%
  as.data.frame() %>%
  setNames(c("f1", "f2", "chao", "pseudo", "turing_p", "turing_b", "n"))
agg <- res %>%
  group_by(n) %>%
  summarise(
    chao = mean(chao),
    pseudo = mean(pseudo),
    turing_p = mean(turing_p),
    turing_b = mean(turing_b)
  )

ggplot(data = res, mapping = aes(group = 2 * n / POP_SIZE, y = turing_b)) +
  geom_boxplot() +
  geom_hline(yintercept = POP_SIZE, linetype = "dashed")
