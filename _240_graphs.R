library(igraph)

set.seed(123)
NODES = 1000
PA_EDGES_PER_NODE = c(1, 2, 3)
DRAWS = 20
TRIALS = 5000
BA_FOLDER = "./_200_input/graphs/ba_%s/"

vector_4cliques = function(x) {
  unlist(lapply(cliques(x, min = 4, max = 4), function(y){ y$name }))
}

write_sample = function(s, filename) {
  clean = Filter(length, s)
  if (length(clean) > 0) {
    writeLines(
      unlist(lapply(clean, paste, collapse=",")),
      con = filename
    )
  }
}

produce_samples = function(
    graph,
    roots,
    order,
    filenames
    ){
  make_ego_graph(graph, order = order, nodes = roots) %>%
  lapply(function(x) { V(x)$name }) %>%
  write_sample(filenames[1])
}

# compute_block_matrix = function(block_nodes, block_edges){
#   p11 = block_edges[1] * 2 / (block_nodes[1] * (block_nodes[1] - 1))
#   p22 = block_edges[2] * 2 / (block_nodes[2] * (block_nodes[2] - 1))
#   p12 = block_edges[3] / (block_nodes[1] * block_nodes[2])
#   rbind(c(p11, p12), c(p12, p22))
# }

for (e in PA_EDGES_PER_NODE) {
  G =
    sample_pa(NODES, 1, m = e, directed = FALSE) %>%
    set.vertex.attribute("name", value = seq(1, NODES))
  c(
    "N",
    "edges",
    "triangles",
    "4_cliques"
    ) %>%
  write(
    file = sprintf(paste(BA_FOLDER, "metadata.csv", sep = ""), e),
    ncolumns = 4,
    append = TRUE,
    sep = ","
    )
  c(
    vcount(G),
    ecount(G),
    sum(count_triangles(G)) / 3,
    length(cliques(G, min = 4, max = 4))
    ) %>%
  write(
    file = sprintf(paste(BA_FOLDER, "metadata.csv", sep = ""), e),
    ncolumns = 4,
    append = TRUE,
    sep = ","
    )
  for (t in seq(1, TRIALS)) {
    roots = sample(1:NODES, DRAWS)
    filenames = c(
      sprintf(paste(BA_FOLDER, "nodes_%s.csv", sep = ""), e, t)
      )
    produce_samples(G, roots, 2, filenames)
  }
}
