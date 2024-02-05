import random

from joblib import Parallel, delayed

import igraph as ig 
import numpy as np

N = 10000
DENSITIES = [0.05, 0.1, 0.2]
R = 1
Q = 0.25
TRIALS = 1000
DRAWS = 50
SEED = 777
DATA_FOLDER = "./_200_input/graphs/ba_{}/graph_{}.csv"


def degbias_ffs2(graph, q, r, root_name, rng):
    """
    Perform a degree-biased forest-fire 2 wave sampling starting from a given
    root node.

    This function implements a degree-biased forest-firs sampling algorithm, where
    nodes are selected based on their degrees, with higher-degree nodes being
    more likely to be chosen. The sampling starts from a specified root node
    and proceeds in 2 waves, with each wave selecting a number of neighbors
    proportional to their degrees.

    Args:
        graph (igraph.Graph): The graph object on which the sampling is performed.

        q (float): The probability parameter for the negative binomial distribution
        used to determine the number of neighbors to select in each wave.

        r (int): The number of trials parameter for the negative binomial distribution.

        root_name (str): The identifier of the root node from which the sampling begins.

        rng (numpy.random.Generator): A random number generator instance for
        reproducible results.

    Returns:
        set: A set of vertices representing the sample obtained through the
        degree-biased forest-fire sampling process.

    Example:
        >>> import igraph as ig
        >>> import numpy as np
        >>> rng = np.random.default_rng(seed=777)
        >>> G = ig.Graph.Erdos_Renyi(10, 0.5, directed=False)
        >>> G.vs["name"] = [str(i) for i in range(1, G.vcount() + 1)]
        >>> sample = degbias_ffs2(G, 0.25, 1, "1", rng)
        >>> print(sample)
        {...} # Set of vertices included in the sample
    """
    # Find the root vertex in the graph
    root = graph.vs.find(root_name)

    # Get the neighbors of the root vertex
    first_ngbhd = root.neighbors()

    # Determine the number of neighbors to select in the first wave
    # Calculate the degrees of the neighbors and sample based on degrees
    n = min([rng.negative_binomial(r, q) + 1, len(first_ngbhd)])
    first_degs = [graph.vs.find(i["name"]).degree() for i in first_ngbhd]
    first_wave = set(rng.choice(
        first_ngbhd,
        size = n,
        replace = False,
        p = np.array(first_degs) / sum(first_degs),
        shuffle = False
    ))
    # Initialize the sample with the first wave of neighbors and the root
    first_wave.add(root)
    sample = set(first_wave)

    # Apply degree-based sampling for each of neighbours of nodes from first wave
    for v in first_wave:
        # Get the neighbors of the current vertex excluding those already in the first wave!
        nghbd = set(v.neighbors()).difference(first_wave)
        if len(nghbd) == 0 or v == root:
            continue
        degs = [graph.vs.find(i["name"]).degree() for i in nghbd]
        n = min([rng.negative_binomial(r, q) + 1, len(nghbd)])
        wave = rng.choice(
            list(nghbd),
            size = n,
            replace = False,
            p = np.array(degs) / sum(degs),
            shuffle = False
        )
        sample.update(wave)
    return sorted(list(sample))


def write_samples(samples, filename):
    """
    Write a collection of samples to a file, with each sample on a separate line.

    Each sample is a list of strings, and each string within a sample represents
    a vertex in the graph. The samples are written to the file in the order they
    appear in the input list, with each sample separated by a newline character.
    Within each sample, the vertices are comma-separated.

    Args:
        samples (list): A list of lists, where each inner list represents a sample
        of vertices.

        filename (str): The path to the file where the samples will be written.

    Example:
        >>> samples = [['A', 'B'], ['C', 'D', 'E']]
        >>> filename = 'output.txt'
        >>> write_samples(samples, filename)
        # Writes to 'output.txt':
        # A,B
        # C,D,E
    """
    with open(filename, "a") as f:
        for s in samples:
            for i in s:
                f.write(i)
                if i == s[-1]:
                    continue
                f.write(',')
            f.write('\n')


def simulate_trial(h, t, seed):
    Ne = 0.5 * h * N * (N - 1)
    rng = np.random.default_rng(seed=seed)
    random.seed(seed)
    G = ig.Graph().Barabasi(N, m = round(Ne / N), directed = False)
    G.vs["name"] = [str(i) for i in range(1, G.vcount() + 1)]
    samples = []
    root_names = rng.choice(
        range(1, N + 1),
        p = np.array(G.degree()) / sum(G.degree()),
        size = DRAWS,
        replace = False,
        shuffle = False
    )
    for name in root_names:
        S = degbias_ffs2(G, Q, R, str(name), rng)
        samples.append([i["name"] for i in S])
    dataset_name = DATA_FOLDER.format(h, t + 1)
    write_samples(samples, dataset_name)
    return None


def main():
    for h in DENSITIES:
        Parallel(n_jobs=-1)(delayed(simulate_trial)(h, t, SEED) for t in range(TRIALS))


if __name__ == "__main__":
    main()
