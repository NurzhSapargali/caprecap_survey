import random

import igraph as ig
import pandas as pd
import numpy as np


def random_walk_sampling_wor_weighted(G, start_node, max_nodes):
    nodes = set([start_node])
    current_node = start_node
    while len(nodes) < max_nodes:
        #print(current_node["name"])
        neighbors = set([G.vs[i] for i in G.neighbors(current_node, mode = "out")])
        #neighbors = list(neighbors - nodes)
        neighbors = list(neighbors)
        if len(neighbors) == 0:
            break
        current_name = current_node["name"]
        weights = [G.es[G.get_eid(current_name, i["name"])]["w"] for i in neighbors]
        next_node = random.choices(neighbors, weights = weights, k = 1)[0]
        nodes.add(next_node)
        current_node = next_node
    return nodes

adj_df = pd.read_csv("./_900_output/data/user_nets/user_net_254.csv")
adj_df = adj_df[adj_df["i"] != adj_df["j"]]
adj_df[['i', 'j']] = adj_df[['i', 'j']].astype(str)
G = ig.Graph.DataFrame(
    adj_df,
    directed = True,
    use_vids = False
)
G.vs["name"] = [str(i) for i in range(G.vcount())]

rng = np.random.default_rng(111)
roots = rng.choice(G.vcount(), replace = False, size = 10000)
for i in roots:
    #S = set(G.random_walk(i, 8000, mode = "out"))
    S = G.neighborhood(vertices = i, order = 10, mode = "out")
    S = [G.vs[i]["name"] for i in S]
    with open("./_200_input/graphs/user_net_254/graph_1.csv", "a") as f:
        for i,s in enumerate(S):
            f.write(s)
            if i != len(S) - 1:
                f.write(",")
            else:
                f.write("\n")
    
yes = [list(random_walk_sampling_wor_weighted(G, G.vs[i], 8000)) for i in range(100000)]