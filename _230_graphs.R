library(igraph)

set.seed(123)
NODES <- 1000
PA_EDGES_PER_NODE <- c(1, 2, 3)
ER_EDGES <- c(1000, 2000, 3000)
DRAWS <- 20
TRIALS <- 1000
BA_FOLDER <- "./_200_input/graphs/ba_%s/"
ER_FOLDER <- "./_200_input/graphs/er_%s/"

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

produce_samples <- function(
    graph,
    roots,
    order,
    filenames
    ){
  sample_graphs <- make_ego_graph(graph, order = order, nodes = roots)
  
  sample_graphs %>% 
    lapply(function(x) { V(x)$name }) %>%
    write_sample(filenames[1])
       # paste(DATA_FOLDER, "nodes_%s.csv", sep = ""), e, t)
  # 
  # sample_graphs %>%
  #   lapply(as_edgelist) %>%
  #   lapply(t) %>%
  #   lapply(as.vector) %>%
  #   write_sample(filenames[2])
  #       # paste(DATA_FOLDER, "edges_%s.csv", sep = ""), e, t)
  # 
  # sample_graphs %>%
  #   lapply(triangles) %>%
  #   lapply(function(x){ x$name })  %>%
  #   write_sample(filenames[3])
  #       # paste(DATA_FOLDER, "tris_%s.csv", sep = ""), e, t)
  # 
  # sample_graphs %>%
  #   lapply(vector_4cliques) %>%
  #   write_sample(filenames[4])
  #     # paste(DATA_FOLDER, "4cliques_%s.csv", sep = ""), e, t)
}

# compute_block_matrix <- function(block_nodes, block_edges){
#   p11 <- block_edges[1] * 2 / (block_nodes[1] * (block_nodes[1] - 1))
#   p22 <- block_edges[2] * 2 / (block_nodes[2] * (block_nodes[2] - 1))
#   p12 <- block_edges[3] / (block_nodes[1] * block_nodes[2])
#   rbind(c(p11, p12), c(p12, p22))
# }

for (t in seq(1, TRIALS)) {
  roots <- sample(1:NODES, DRAWS)
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
    write(
      meta,
      file = sprintf(paste(BA_FOLDER, "metadata.csv", sep = ""), e),
      append = TRUE
      )
    filenames <- c(
      sprintf(paste(BA_FOLDER, "nodes_%s.csv", sep = ""), e, t)
      # sprintf(paste(BA_FOLDER, "edges_%s.csv", sep = ""), e, t),
      # sprintf(paste(BA_FOLDER, "tris_%s.csv", sep = ""), e, t),
      # sprintf(paste(BA_FOLDER, "4cliques_%s.csv", sep = ""), e, t)
    )
    produce_samples(G, roots, 2, filenames)
  }
  # for (e in ER_EDGES){
  #   G <-
  #     sample_gnm(NODES, e, directed = FALSE) %>%
  #     set.vertex.attribute("name", value = seq(1, NODES))
  #   meta <- c(
  #     vcount(G),
  #     ecount(G),
  #     sum(count_triangles(G)) / 3,
  #     length(cliques(G, min = 4, max = 4))
  #   )
  #   write(
  #     meta,
  #     file = sprintf(paste(ER_FOLDER, "metadata.csv", sep = ""), e),
  #     append=TRUE
  #     )
  #   filenames <- c(
  #     sprintf(paste(ER_FOLDER, "nodes_%s.csv", sep = ""), e, t),
  #     sprintf(paste(ER_FOLDER, "edges_%s.csv", sep = ""), e, t),
  #     sprintf(paste(ER_FOLDER, "tris_%s.csv", sep = ""), e, t),
  #     sprintf(paste(ER_FOLDER, "4cliques_%s.csv", sep = ""), e, t)
  #   )
  #   produce_samples(G, roots, 2, filenames)
  #   for (i in seq_along(SBM_BLOCK_EDGE_PROPS)){
  #     K <- compute_block_matrix(SBM_BLOCK_SIZES, SBM_BLOCK_EDGE_PROPS[[i]] * e)
  #     G <- 
  #       sample_sbm(NODES, K, SBM_BLOCK_SIZES, directed = FALSE, loops = FALSE) %>%
  #       set.vertex.attribute("name", value = seq(1, NODES))
  #     meta <- c(
  #       vcount(G),
  #       ecount(G),
  #       sum(count_triangles(G)) / 3,
  #       length(cliques(G, min = 4, max = 4))
  #     )
  #     gtype <- names(SBM_BLOCK_EDGE_PROPS)[i]
  #     write(
  #       meta,
  #       file = sprintf(paste(SBM_FOLDER, "metadata.csv", sep = ""), e, gtype),
  #       append=TRUE
  #       )
  #     filenames <- c(
  #       sprintf(paste(SBM_FOLDER, "nodes_%s.csv", sep = ""), e, gtype, t),
  #       sprintf(paste(SBM_FOLDER, "edges_%s.csv", sep = ""), e, gtype, t),
  #       sprintf(paste(SBM_FOLDER, "tris_%s.csv", sep = ""), e, gtype, t),
  #       sprintf(paste(SBM_FOLDER, "4cliques_%s.csv", sep = ""), e, gtype, t)
  #     )
  #     produce_samples(G, roots, 2, filenames)
  #   }
  # }
}
