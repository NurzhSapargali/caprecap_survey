library(igraph)

set.seed(123)
NODES = c(5000, 10000)
DRAWS = 50
TRIALS = 1
BA_FOLDER = "./_200_input/graphs/ba_%s/"

# Write a sample graph to a file
# Inputs:
#  s: The graph sample to write
#  filename: The name of the file to write to
# Outputs: None
# Side Effects: Writes to a file
write_sample <- function(s, filename) {
  # Remove empty elements from the sample
  clean <- Filter(length, s)
  # Only write to the file if there are elements in the cleaned sample
  if (length(clean) > 0) {
    # Write the cleaned sample to the file
    writeLines(
      unlist(lapply(clean, paste, collapse=",")),
      con = filename
    )
  }
}

# Produce snowball samples from a given graph
# Inputs:
# graph: The graph to produce samples from
# roots: The starting vertices from which to snowball
# order: The neighbourhood radius of snowball sample
# filenames: The filenames to write the samples to
# Outputs: None
# Side Effects: Writes to a file
produce_samples <- function(graph, roots, order, filenames) {
  # Create an ego graph from the input graph
  ego_graph <- make_ego_graph(graph, order = order, nodes = roots)
  # Extract the vertex names from the ego graph
  vertex_names <- lapply(ego_graph, function(x) { V(x)$name })
  # Write the vertex names to the specified file
  write_sample(vertex_names, filenames[1])
}

sbm_block_sizes = lapply(NODES, function(x) {c(x * 0.5, x * 0.35, x * 0.15)})
pa_edges_per_node = lapply(
  NODES,
  function(x) {
    round(c(
      (x - 1) / 2 * 0.20,
      (x - 1) / 2 * 0.40
    ))
  }
)
names(sbm_block_sizes) = NODES
names(pa_edges_per_node) = NODES
sbm_matrices = list(
  matrix(
    c(0.15, 0.05, 0.05, 0.05, 0.15, 0.05, 0.05, 0.05, 0.15),
    nrow = 3,
    ncol = 3,
    byrow = TRUE
  ),
  matrix(
    c(0.6, 0.08, 0.08, 0.08, 0.6, 0.08, 0.08, 0.08, 0.6),
    nrow = 3,
    ncol = 3,
    byrow = TRUE
  )
)
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
