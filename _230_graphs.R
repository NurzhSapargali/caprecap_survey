library(igraph)

# Simulation constants
set.seed(123)
NODES = 10000
EDGES = c(70000, 100000)
SBM_BLOCK_SIZES = c(0.5 * NODES, 0.3 * NODES, 0.2 * NODES)
# prob_mat = matrix(
#   c(
#     0.004345726,
#     0.002172863,
#     0.002172863,
#     0.002172863,
#     0.004345726,
#     0.002172863,
#     0.002172863,
#     0.002172863,
#     0.004345726
#   ),
#   byrow = TRUE,
#   nrow = 3,
#   ncol = 3
# )
# prob_mat = matrix(
#   c(
#     0.002897152,
#     0.001448576,
#     0.001448576,
#     0.001448576,
#     0.002897152,
#     0.001448576,
#     0.001448576,
#     0.001448576,
#     0.002897152
#   ),
#   byrow = TRUE,
#   nrow = 3,
#   ncol = 3
# )
SBM_BLOCK_MATRICES = list(
  matrix(
    c(0.002, 0.001, 0.001, 0.001, 0.002, 0.001, 0.001, 0.001, 0.002),
    byrow = TRUE,
    nrow = 3,
    ncol = 3
  ),
  matrix(
    c(0.0027, 0.0014, 0.0014, 0.0014, 0.0027, 0.0014, 0.0014, 0.0014, 0.0027),
    byrow = TRUE,
    nrow = 3,
    ncol = 3
  )
)
DRAWS = 50
TRIALS = 500
BA_FOLDER = "./_200_input/graphs/ba_%s/"
SBM_FOLDER = "./_200_input/graphs/sbm_%s/"

# Write a sample graph to a file
# Inputs:
#  s: The graph sample to write
#  filename: The name of the file to write to
# Outputs: None
# Side Effects: Writes to a file
write_sample = function(s, filename) {
  # Remove empty elements from the sample
  clean = Filter(length, s)
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
produce_samples = function(graph, roots, order, filenames) {
  # Create an ego graph from the input graph
  ego_graph = make_ego_graph(graph, order = order, nodes = roots)
  # Extract vertex names from the ego graph
  vertex_names = lapply(ego_graph, function(x) { V(x)$name })
  # Write vertex names to the specified file
  write_sample(vertex_names, filenames[1])
}

# Simulate preferential attachment graphs
for (e in EDGES) {
  write(
    c("N", "edges"),
    file = sprintf(
      paste(BA_FOLDER, "metadata_%s.csv", sep = ""),
      e / 1000,
      e / 1000
    ),
    ncolumns = 2,
    append = TRUE,
    sep = ","
  )
  for (t in seq(1, TRIALS)) {
    gpa =
      sample_pa(NODES, 1, m = e / NODES, directed = FALSE) %>%
      set.vertex.attribute("name", value = seq(1, NODES))
    # Write out number of edges and nodes of simulated graph
    write(
      c(vcount(gpa), ecount(gpa)),
      file = sprintf(
        paste(BA_FOLDER, "metadata_%s.csv", sep = ""),
        e / 1000,
        e / 1000
      ),
      ncolumns = 2,
      append = TRUE,
      sep = ","
    )
    # Write out samples of nodes from snowball sampling to a separate file
    roots = sample(1:NODES, DRAWS)
    filenames = c(
      sprintf(paste(BA_FOLDER, "graph_%s.csv", sep = ""), e / 1000, t)
    )
    produce_samples(gpa, roots, 1, filenames)
  }
}
# Simulate stochastic block model graphs
for (m in seq_along(SBM_BLOCK_MATRICES)) {
  write(
    c("N", "edges"),
    file = sprintf(
      paste(SBM_FOLDER, "metadata_%s.csv", sep = ""),
      EDGES[m] / 1000,
      EDGES[m] / 1000
    ),
    ncolumns = 2,
    append = TRUE,
    sep = ","
  )
  for (t in seq(1, TRIALS)) {
    gsbm =
      sample_sbm(
        NODES,
        SBM_BLOCK_MATRICES[[m]],
        SBM_BLOCK_SIZES,
        directed = FALSE
      ) %>%
      set.vertex.attribute("name", value = seq(1, NODES))
    # Write out number of edges and nodes of simulated graph
    write(
      c(vcount(gsbm), ecount(gsbm)),
      file = sprintf(
        paste(SBM_FOLDER, "metadata_%s.csv", sep = ""),
        EDGES[m] / 1000,
        EDGES[m] / 1000
      ),
      ncolumns = 2,
      append = TRUE,
      sep = ","
    )
    # Write out samples of nodes from snowball sampling to a separate file
    roots = sample(1:NODES, DRAWS)
    filenames = c(
      sprintf(paste(SBM_FOLDER, "graph_%s.csv", sep = ""), EDGES[m] / 1000, t)
    )
    produce_samples(gsbm, roots, 1, filenames)
  }
}

