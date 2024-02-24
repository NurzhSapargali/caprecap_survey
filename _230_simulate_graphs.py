#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
File: 230_simulate_graphs.py
---------------------------
This script simulates a series of graph sampling experiments. It generates
Barabasi-Albert graphs with specified densities and applies a degree-biased
forest-fire sampling algorithm to sample vertices from the graph. The results
of the sampling process are written to a file.

Author: Nurzhan Sapargali
Date created: 2024-02-24
"""

import random

from joblib import Parallel, delayed

import igraph as ig 
import numpy as np

N = 10000 # True population
DENSITIES = [0.05, 0.1, 0.2] # Graph densities
R = 1 # Number of trials for negative binomial distribution
Q = 0.25 # Success probability for negative binomial distribution 
TRIALS = 1000 # Number of simulations
DRAWS = 50 # Number of sample draws for each simulation
DATA_FOLDER = "./_200_input/graphs/ba_{}/graph_{}.csv" # Output file name template for samples
METADATA_FILE = "./_200_input/graphs/ba_{}/metadata.csv" # Output file name for metadata


def pareto_sample(p, n, rng):
    """
    Generate a probability-proportional-to-size (PPS) sample given sizes of
    individuals and desirable sample size using the Pareto sampling algorithm.

    Args:
        p (numpy.array): The sizes of each individual for PPS-sampling 

        n (int): The sample size.

        rng (np.random.default_rng): A random number generator for the sampling
        process.

    Returns:
        numpy.array: An array of indices representing the Pareto sample.

    Example:
        >>> p = numpy.array([0.1,  0.2,  0.3,  0.4])
        >>> n =  2
        >>> rng = numpy.random.default_rng(seed=42)
        >>> sample = pareto_sample(p, n, rng)
        >>> print(sample)
        [3 1]
    """
    # Calculate the inclusion probabilities
    pr = p * n
    N = len(pr)
    
    # Generate uniform random numbers
    u = rng.uniform(0,  1, N)
    
    # Calculate the Pareto weights
    q = (u / (1.0 - u)) / (pr / (1.0 - pr))
    
    # Return n indices with smallest weights
    return np.argsort(q)[:n]


def ar_pareto_sample(p, n, rng):
    """
    Generate a Sampford-Rao PPS sample given sizes of individuals and desirable
    sample size using acceptance-rejection Pareto sampling algorithm from:

    Bondesson, L., Imbi Traat, & Anders Lundqvist. (2006). Pareto Sampling versus
    Sampford and Conditional Poisson Sampling. Scandinavian Journal of Statistics,
    33(4), 699–720. http://www.jstor.org/stable/4616953.

    Args:
        p (numpy.array): The sizes of each individual for PPS-sampling 

        n (int): The sample size.

        rng (numpy.random.default_rng): A random number generator for the sampling
        process.

    Returns:
        numpy.array: An array of indices representing a Sampford-Rao sample.

    Example:
        >>> p = numpy.array([0.1,   0.2,   0.3,   0.4])
        >>> n =   2
        >>> rng = numpy.random.default_rng(seed=42)
        >>> sample = ar_pareto_sample(p, n, rng)
        >>> print(sample)
        [3 1]
    """
    # Calculate the inclusion probabilities
    pr = n * p
    N = len(pr)

    # Calculate quantities as described by Bondesson et al. (2006)
    d = sum(pr * (1.0 - pr))
    sigma2_k =  1.0 / (d + pr * (1.0 - pr))
    ckc0 = (1.0 - pr) * sigma2_k**0.5 * np.exp(sigma2_k * pr**2 /  2)
    ckc0 = ckc0 / sum((1.0 - pr) * sigma2_k**0.5 * np.exp(sigma2_k * pr**2 /  2))
    ckc0 = (N - n) * ckc0
    Jk = ckc0 / (1.0 - pr)
    A = min(Jk)

    # Acceptance-rejection filter
    accept = False
    while not accept:
        # Generate a Pareto sample
        S = pareto_sample(p, n, rng)
        pr_s = np.take(pr, S)
        Jk_s = np.take(Jk, S)
        # Calculate the acceptance ratio
        car = sum(1.0 - pr_s) / sum(Jk_s / A * (1.0 - pr_s))
        u = rng.uniform(0,  1)
        # Accept or reject the sample based on the acceptance ratio
        if (u <= car):
            accept = True

    return S


def degbias_ffs2(graph, q, r, root_name, rng):
    """
    Perform a degree-biased forward-first sampling (FFS) starting from a root
    vertex.

    This function implements a degree-biased forward-first sampling algorithm,
    where the sampling process starts from a specified root vertex and proceeds
    by selecting neighbors based on their degrees, with higher-degree vertices
    being more likely to be selected. The algorithm is applied in a two-wave
    manner: in the first wave, neighbors of the root are selected based on their
    degrees, and in the second wave, neighbors of the selected nodes from the
    first wave are selected in the same manner.

    Args:
        graph (igraph.Graph): The graph on which the sampling is performed.

        q (float): The parameter for the negative binomial distribution used to
        determine the number of neighbors to select in each wave.

        r (float): The parameter for the negative binomial distribution used to
        determine the number of neighbors to select in each wave.

        root_name (str): The name of the root vertex from which the sampling
        starts.

        rng (numpy.random.default_rng): A random number generator for the
        sampling process.

    Returns:
        list: A list of vertex names that were selected during the sampling process,
              including the root vertex.

    Example:
        >>> graph = igraph.Graph().Barabasi(100, m=5)
        >>> graph.vs["name"] = [str(i) for i in range(1, graph.vcount() +  1)]
        >>> q =  0.5
        >>> r =  0.5
        >>> root_name = '1'
        >>> rng = numpy.random.default_rng(seed=42)
        >>> selected_vertices = degbias_ffs2(graph, q, r, root_name, rng)
        >>> print(selected_vertices)
        ['1', '10', '15', '2', '24', '26', '29', '34', '41', '5', '6']
    """
    # Find the root vertex in the graph
    root = graph.vs.find(root_name)

    # Get the neighbors of the root vertex
    first_ngbhd = root.neighbors()

    # Determine the number of neighbors to select in the first wave
    # Calculate the degrees of the neighbors and sample based on degrees
    n = min([rng.negative_binomial(r, q) + 1, len(first_ngbhd)])
    first_degs = [graph.vs.find(i["name"]).degree() for i in first_ngbhd]
    p = np.array(first_degs) / sum(first_degs)
    first_wave = np.take(first_ngbhd, ar_pareto_sample(p, n, rng))

    # Initialize the sample with the first wave of neighbors and the root
    sample = set(first_wave).union({root})

    # Apply degree-based sampling for each of neighbours of nodes from first wave
    for v in first_wave:
        # Get the neighbors of the current vertex excluding those already in the first wave!
        ngbhd = set(v.neighbors()).difference(set(first_wave).union({root}))
        ngbhd = sorted(ngbhd)
        if len(ngbhd) == 0:
            continue

        # Produce Sampford-Rao sample on the remaining neighbours based on their degrees
        degs = [graph.vs.find(i["name"]).degree() for i in ngbhd]
        n = min([rng.negative_binomial(r, q) + 1, len(ngbhd)])
        p = np.array(degs) / sum(degs)
        wave = np.take(ngbhd, ar_pareto_sample(p, n, rng))
        sample.update(set(wave))

    return sorted([i["name"] for i in sample])


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


def simulate_trial(h, t):
    """
    Simulate a single trial of a graph sampling experiment.

    This function simulates a single trial of a graph sampling experiment, where
    a Barabasi-Albert graph is generated with a specified density, and a degree-
    biased forest-fire sampling algorithm is applied to sample vertices from
    the graph. The sampling process starts from a set of root vertices, which are
    sampled via PPS sampling based on inclusion probabilities proportional to
    normalised vertex degrees. The results of the sampling process are written
    to a file.

    Args:
        h (float): The density of the Barabasi-Albert graph to be generated.

        t (int): The trial number, used for naming the output file and seeding
        the random number generator.

    Returns:
        None
    """
    # Calculate the number of edges for the Barabasi-Albert graph
    Ne =  0.5 * h * N * (N -  1)

    # Set seed for the random number generator
    random.seed(t + 1)
    rng = np.random.default_rng(t + 1)
    
    # Generate the Barabasi-Albert graph and write it out some metadata
    G = ig.Graph().Barabasi(N, m = round(Ne / N), directed = False)
    G.vs["name"] = [str(i) for i in range(1, G.vcount() +  1)]
    # Unique identifier for the graph as a sanity check of reproducibility
    G["id"] = random.uniform(0, 1)

    # Write the metadata to the file
    meta = [
        G["id"],
        G.vcount(),
        G.density(),
        np.mean(G.degree())
    ]
    write_samples([[str(i) for i in meta]], METADATA_FILE.format(h))
    
    # Initialize the list to store the samples
    samples = []
    
    # Calculate the probability distribution based on the degrees of the vertices
    p = np.array(G.degree()) / sum(G.degree())
    
    # Sample the root vertices
    root_names = np.take(G.vs, ar_pareto_sample(p, DRAWS, rng))
    
    # Perform the degree-biased forest-fire sampling for each root vertex
    for name in [i["name"] for i in root_names]:
        S = degbias_ffs2(G, Q, R, name, rng)
        samples.append(S)
    
    # Format the dataset name
    dataset_name = DATA_FOLDER.format(h, t +  1)
    
    # Write the samples to the file
    write_samples(samples, dataset_name)
    
    return None


def main():
    for h in DENSITIES:
        cols = [
            "id",
            "N",
            "density",
            "avg_degree"
        ]
        write_samples([cols], METADATA_FILE.format(h))
        Parallel(n_jobs=-1)(delayed(simulate_trial)(h, t) for t in range(TRIALS))


if __name__ == "__main__":
    main()
