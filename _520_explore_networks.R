library(igraph)
library(tidyverse)
library(ggpubr)

USER_NETS <- "./_900_output/data/user_nets/"
DATE_IDS <- list(
  "254" = "27 Nov - 28 Nov 2020",
  "255" = "28 Nov - 29 Nov 2020",
  "256" = "29 Nov - 30 Nov 2020",
  "984" = "27 Nov - 28 Nov 2022",
  "985" = "28 Nov - 29 Nov 2022",
  "986" = "29 Nov - 30 Nov 2022"
)

retrieve_date <- function(filename, id_map) {
  for (id in names(id_map)) {
    if (grepl(id, filename, fixed = TRUE)) {
      return(id_map[[id]])
    }
  }
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
filenames <- list.files(USER_NETS, full.names = TRUE)
degs <- list()
for (i in seq_along(filenames)) {
  f <- filenames[i]
  adj_df <- read.csv(f, header = TRUE) %>% filter(i != j)
  g <- graph_from_data_frame(adj_df)
  graph_attr(g, "name") <- retrieve_date(f, DATE_IDS)
  both_degs <- rbind(degree_dist(g, mode = "in"), degree_dist(g, mode = "out"))
  degs[[i]] <- both_degs
}

degs <- do.call(rbind, degs)

ggplot(degs, aes(x = log_degree, y = log_count, color = type)) +
  geom_point() +
  theme_minimal() +
  theme_pubr(base_size = 17) +
  theme(legend.title = element_blank()) +
  xlab("Log degree") +
  ylab("Log count") +
  facet_wrap(
    vars(date),
    ncol = 3,
    nrow = 2
  )
ggsave(
  "test.pdf",
  width = 297,
  height = 165,
  units = "mm"
)