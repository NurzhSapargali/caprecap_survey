library(igraph)

set.seed(123)
NODES <- 1000
PA_EDGES_PER_NODE <- c(1, 2, 3)
DRAWS <- 20
DATA_FOLDER <- "./_200_input/graphs/ba_%s/"
TRIALS <- 25

vector_4cliques <- function(x) {
  unlist(lapply(cliques(x, min = 4, max = 4), function(y){ y$name }))
}

write_sample <- function(s, filename) {
  clean <- Filter(length, s)
  if (length(clean) > 0) {
    writeLines(
      unlist(lapply(clean, paste, collapse=",")),
      con = filename
    )
  }
}

for (e in PA_EDGES_PER_NODE) {
  G <-
    sample_pa(NODES, 1, m = e, directed = FALSE) %>%
    set.vertex.attribute("name", value = seq(1, NODES))
  meta <- c(
    vcount(G),
    ecount(G),
    sum(count_triangles(G)) / 3,
    length(cliques(G, min = 4, max = 4))
    )
  write(meta, file = sprintf(paste(DATA_FOLDER, "metadata.csv", sep = ""), e))
  for (t in seq(1, TRIALS)) {
    roots <- sample(1:NODES, DRAWS)
    sample_graphs <- make_ego_graph(G, order = 2, nodes = roots)

    sample_graphs %>%
      lapply(function(x) { V(x)$name }) %>%
      write_sample(
        sprintf(paste(DATA_FOLDER, "nodes_%s.csv", sep = ""), e, t)
      )

    sample_graphs %>%
      lapply(as_edgelist) %>%
      lapply(t) %>%
      lapply(as.vector) %>%
      write_sample(
        sprintf(paste(DATA_FOLDER, "edges_%s.csv", sep = ""), e, t)
      )

    sample_graphs %>%
      lapply(triangles) %>%
      lapply(function(x){ x$name })  %>%
      write_sample(
        sprintf(paste(DATA_FOLDER, "tris_%s.csv", sep = ""), e, t)
      )

    sample_graphs %>%
      lapply(vector_4cliques) %>%
      write_sample(
        sprintf(paste(DATA_FOLDER, "4cliques_%s.csv", sep = ""), e, t)
      )
  }
}
