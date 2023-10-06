library(igraph)

set.seed(123)
NODES = c(1000, 5000, 10000)
PA_EDGES_PER_NODE = c(1, 2, 3)
DRAWS = 20
TRIALS = 500
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

for (n in NODES) {
  for (e in PA_EDGES_PER_NODE) {
    G =
      sample_pa(n, 1, m = e, directed = FALSE) %>%
      set.vertex.attribute("name", value = seq(1, n))
    c(
      "N",
      "edges",
      "triangles",
      "4_cliques"
      ) %>%
    write(
      file = sprintf(paste(BA_FOLDER, "metadata_%s.csv", sep = ""), e, n),
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
      file = sprintf(paste(BA_FOLDER, "metadata_%s.csv", sep = ""), e, n),
      ncolumns = 4,
      append = TRUE,
      sep = ","
      )
    for (t in seq(1, TRIALS)) {
      roots = sample(1:n, DRAWS)
      filenames = c(
        sprintf(paste(BA_FOLDER, "nodes_%s_%s.csv", sep = ""), e, t, n)
        )
      produce_samples(G, roots, 2, filenames)
    }
  }
}
