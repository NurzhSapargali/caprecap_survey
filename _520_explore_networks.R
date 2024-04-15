library(igraph)
library(dplyr)
library(ggplot2)
library(ggpubr)

USER_NETS <- "./_900_output/data/user_nets/"
FIGURE_FOLDER <- "./_900_output/figures/user_nets/"
DATE_IDS <- list(
  "250" = "23 - 24 Nov 2020",
  "251" = "24 - 25 Nov 2020",
  "252" = "25 - 26 Nov 2020",
  "253" = "26 - 27 Nov 2020",
  "254" = "27 - 28 Nov 2020",
  "980" = "23 - 24 Nov 2022",
  "981" = "24 - 25 Nov 2022",
  "982" = "25 - 26 Nov 2022",
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
  comps <- components(g, mode = "weak")
  w_comps <- comps$no
  avg_deg <- m / n
  rel_size <- max(comps$csize) / n
  return(c(n, m, avg_deg, clust_coef, w_comps, rel_size))
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
nov_2020 <- list()
nov_2022 <- list()
graph_summaries <- list()
for (i in seq_along(filenames)) {
  f <- filenames[i]
  date_range <- retrieve_date(f, DATE_IDS)
  if (is.null(date_range)) {
    next
  }

  adj_df <- read.csv(f, header = TRUE) %>% dplyr::filter(i != j)
  g <- graph_from_data_frame(adj_df)

  graph_summaries[[date_range]] <- summarize_graph(g)
  if (grepl("2020", date_range)) {
    nov_2020[[date_range]] <- adj_df
  } else {
    nov_2022[[date_range]] <- adj_df
  }
}

graph_summaries <- do.call(rbind, graph_summaries)
graph_summaries <- as.data.frame(graph_summaries)
colnames(graph_summaries) <- c(
  "vertices",
  "edges",
  "avg_deg",
  "clust_coef",
  "w_comps",
  "rel_size"
)
print(graph_summaries)

nov_2020 <- do.call(rbind, nov_2020) %>%
  distinct(i, j, .keep_all = TRUE) %>%
  graph_from_data_frame()
graph_attr(nov_2020, "name") <- "Nov 2020"

nov_2022 <- do.call(rbind, nov_2022) %>%
  distinct(i, j, .keep_all = TRUE) %>%
  graph_from_data_frame()
graph_attr(nov_2022, "name") <- "Nov 2022"

degs <- list(
  degree_dist(nov_2020, mode = "in"),
  degree_dist(nov_2020, mode = "out"),
  degree_dist(nov_2022, mode = "in"),
  degree_dist(nov_2022, mode = "out")
)
degs <- do.call(rbind, degs)

ggplot(degs, aes(x = log_degree, y = log_count, color = type)) +
  geom_point(size = 1.5, alpha = 0.5) +
  theme_minimal() +
  theme_pubr(base_size = 15) +
  theme(
    legend.title = element_blank(),
    legend.text = element_text(size = 15),
    text = element_text(size = 15),
    strip.text = element_text(size = 15)
  ) +
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) +
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

merged_summary <- rbind(
  summarize_graph(nov_2020),
  summarize_graph(nov_2022)
) %>%
  as.data.frame()
colnames(merged_summary) <- colnames(graph_summaries)
print(merged_summary)
