library(igraph)
library(dplyr)
library(ggplot2)
library(ggpubr)

USER_NETS <- "./_900_output/data/user_nets/"
FIGURE_FOLDER <- "./_900_output/figures/user_nets/"
DATE_IDS <- list(
  "254" = "27 - 28 Nov 2020",
  "255" = "28 - 29 Nov 2020",
  "256" = "29 - 30 Nov 2020",
  "984" = "27 - 28 Nov 2022",
  "985" = "28 - 29 Nov 2022",
  "986" = "29 - 30 Nov 2022"
)

retrieve_date <- function(filename, id_map) {
  for (id in names(id_map)) {
    if (grepl(id, filename, fixed = TRUE)) {
      return(id_map[[id]])
    }
  }
}

summarize_graph <- function(g) {
  n <- vcount(g)
  m <- ecount(g)
  clust_coef <- transitivity(g, type = "global")
  w_comps <- count_components(g, mode = "weak")
  s_comps <- count_components(g, mode = "strong")
  avg_deg <- mean(degree(g, mode = "in"))
  return(c(n, m, avg_deg, clust_coef, w_comps, s_comps))
}

degree_dist <- function(g, mode = "in") {
  deg_dist <- degree(g, mode = mode)
  deg_dist <- table(deg_dist[deg_dist > 0]) %>% as.data.frame()
  colnames(deg_dist) <- c("degree", "count")
  deg_dist$degree <- as.numeric(as.character(deg_dist$degree))
  deg_dist <- mutate(deg_dist, log_count = log(count), log_degree = log(degree))
  deg_dist$type <- paste(mode, "degree", sep = "-")
  deg_dist$date <- as.factor(graph_attr(g, "name"))
  deg_dist
}

# Load the networks
filenames <- list.files(USER_NETS, full.names = TRUE, pattern = "\\.csv$")
degs <- list()
graph_summaries <- list()
weights <- list()
for (i in seq_along(filenames)) {
  f <- filenames[i]
  adj_df <- read.csv(f, header = TRUE) %>% dplyr::filter(i != j)
  g <- graph_from_data_frame(adj_df)

  date_range <- retrieve_date(f, DATE_IDS)

  weight_dist <- table(adj_df$w[adj_df$w > 0]) %>% as.data.frame()
  colnames(weight_dist) <- c("w", "count")
  weight_dist$w <- as.numeric(as.character(weight_dist$w))
  weight_dist <- mutate(weight_dist, log_w = log(w), log_count = log(count))
  weight_dist$date <- as.factor(date_range)
  weights[[i]] <- weight_dist

  graph_attr(g, "name") <- date_range
  both_degs <- rbind(degree_dist(g, mode = "in"), degree_dist(g, mode = "out"))
  degs[[i]] <- both_degs

  graph_summaries[[date_range]] <- summarize_graph(g)
}

degs <- do.call(rbind, degs)

ggplot(degs, aes(x = log_degree, y = log_count, color = type)) +
  geom_point(size = 1.5, alpha = 0.5) +
  theme_minimal() +
  theme_pubr(base_size = 15) +
  theme(legend.title = element_blank()) +
  xlab("Log degree") +
  ylab("Log count") +
  facet_wrap(
    vars(date),
    ncol = 3,
    nrow = 2
  )
ggsave(
  paste(FIGURE_FOLDER, "deg_dist.pdf", sep = ""),
  width = 297,
  height = 165,
  units = "mm"
)

weights <- do.call(rbind, weights)

ggplot(weights, aes(x = log_w, y = log_count)) +
  geom_point(size = 1.5, alpha = 0.5) +
  theme_minimal() +
  theme_pubr(base_size = 15) +
  theme(legend.title = element_blank()) +
  xlab("Log mentions") +
  ylab("Log count") +
  facet_wrap(
    vars(date),
    ncol = 3,
    nrow = 2
  )
ggsave(
  paste(FIGURE_FOLDER, "weight_dist.pdf", sep = ""),
  width = 297,
  height = 165,
  units = "mm"
)

graph_summaries <- do.call(rbind, graph_summaries)
graph_summaries <- as.data.frame(graph_summaries)
colnames(graph_summaries) <- c(
  "vertices",
  "edges",
  "avg_deg",
  "clust_coef",
  "w_comps",
  "s_comps"
)
print(graph_summaries)