#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 15:43:55 2022

@author: nurzhan
"""
from glob import glob

import numpy as np
import pandas as pd


ESTIMATES_DATA = "./_900_output/data/"
ALPHAS = [1.0, 5.0, 10.0, None]
N = 1000

for alpha in ALPHAS:
    stdevs = []
    means = []
    if alpha:
        data = pd.read_csv(ESTIMATES_DATA + "diffp/estimates_{}.csv".format(alpha),
                           names=["alpha_hat", "N_u", "N_o", "T", "alpha"])
        bench = pd.read_csv(ESTIMATES_DATA + "diffp/benchmarks_{}.csv".format(alpha),
                            names=["N_hat", "N_o", "T", "alpha", "estimator"])
    else:
        data = pd.read_csv(ESTIMATES_DATA + "eqp/estimates.csv",
                           names=["alpha_hat", "N_u", "N_o", "T"])
        bench = pd.read_csv(ESTIMATES_DATA + "eqp/benchmarks.csv",
                            names=["N_hat", "N_o", "T", "estimator"])
    data["N_hat"] = data["N_u"] + data["N_o"]
    bench["N_hat"] = pd.to_numeric(bench["N_hat"], errors="coerce")
    data["estimator"] = "Pseudolikelihood"
    df = pd.concat([data, bench])
    for t in df["T"].unique():
        cut_mean = df[df["T"] == t].groupby("estimator").mean()[["N_hat"]]
        cut_mean["T"] = t
        cut_mean["alpha"] = alpha
        cut_var = df[df["T"] == t].groupby("estimator").std()[["N_hat"]]
        cut_var["T"] = t
        cut_var["alpha"] = alpha
        means.append(cut_mean)
        stdevs.append(cut_var)
    means = pd.concat(means)
    stdevs = pd.concat(stdevs)
    if alpha:
        means.to_csv(ESTIMATES_DATA + "diffp/sim_means_{}.csv".format(alpha))
        stdevs.to_csv(ESTIMATES_DATA + "diffp/sim_stdevs_{}.csv".format(alpha))
    else:
        means.to_csv(ESTIMATES_DATA + "eqp/sim_means.csv")
        stdevs.to_csv(ESTIMATES_DATA + "eqp/sim_stdevs.csv")
    