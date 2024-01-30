import igraph as ig 
import numpy as np

N = 10000
DENSITIES = [0.05, 0.1, 0.2]
R = 1
Q = 0.25
TRIALS = 100
DRAWS = 50
SEED = 777
DATA_FOLDER = "./_200_input/graphs/ba_{}/graph_{}.csv"


def write_samples(samples, filename):
    with open(filename, "a") as f:
        for s in samples:
            for i in s:
                f.write(i)
                if i == s[-1]:
                    continue
                f.write(',')
            f.write('\n')


def degbias_ffs2(graph, q, r, root_name, rng):
    root = graph.vs.find(root_name)
    first_ngbhd = root.neighbors()
    n = min([rng.negative_binomial(r, q) + 1, len(first_ngbhd)])
    first_degs = [graph.vs.find(i["name"]).degree() for i in first_ngbhd]
    first_wave = rng.choice(
        first_ngbhd,
        size = n,
        replace = False,
        p = np.array(first_degs) / sum(first_degs)
    )
    sample = set(first_wave)
    sample.add(root)
    for v in first_wave:
        nghbd = set(v.neighbors())
        nghbd = nghbd.difference(first_wave)
        degs = [graph.vs.find(i["name"]).degree() for i in nghbd]
        n = min([rng.negative_binomial(r, q) + 1, len(nghbd)])
        wave = rng.choice(
            list(nghbd),
            size = n,
            replace = False,
            p = np.array(degs) / sum(degs)
        )
        sample.update(wave)
    return sample


def main():
    rng = np.random.default_rng(SEED)
    for h in DENSITIES:
        Ne = 0.5 * h * N * (N - 1)
        for t in range(TRIALS):
            G = ig.Graph().Barabasi(N, m = round(Ne / N), directed = False)
            G.vs["name"] = [str(i) for i in range(1, G.vcount() + 1)]
            samples = []
            root_names = rng.choice(
                range(1, N + 1),
                p = np.array(G.degree()) / sum(G.degree()),
                size = DRAWS,
                replace = False
            )
            for name in root_names:
                S = degbias_ffs2(G, Q, R, str(name), rng)
                samples.append([i["name"] for i in S])
            dataset_name = DATA_FOLDER.format(h, t + 1)
            write_samples(samples, dataset_name)


if __name__ == "__main__":
    main()
