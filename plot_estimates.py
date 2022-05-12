#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 15:43:55 2022

@author: nurzhan
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

chaos = pd.read_csv("chaos.csv", header=None)
chaos.columns = ["N_hat", "alpha", "T"]
chaos["estimator"] = "Chao"
chaos_corr = pd.read_csv("chaos_corr.csv", header=None)
chaos_corr.columns = ["N_hat", "alpha", "T"]
chaos_corr["estimator"] = "Chao (corrected)"
links = pd.read_csv("links.csv", header=None)
links.columns = ["N_hat", "alpha", "T"]
links["estimator"] = "Lincoln & Schnabel"
schnab = pd.read_csv("schnab.csv", header=None)
schnab.columns = ["N_hat", "alpha", "T"]
schnab["estimator"] = "Lincoln & Schnabel"
estis = pd.read_csv("estis.csv", header=None)
estis.columns = ["alpha_hat", "N_hat", "alpha", "T"]
estis["estimator"] = "Pseudo-likelihood"
out = pd.concat([links, schnab, chaos, chaos_corr, estis.iloc[:,1:]])
out["1 / N_hat"] = 1 / out["N_hat"]
for alpha in out["alpha"].unique():
    fig, ax = plt.subplots(figsize=(18,12))
    sns.boxplot(x="1 / N_hat", y="T", hue="estimator", orient="h",
                data=out[out["alpha"] == alpha], whis=[5, 95], meanline=True,
                ax=ax)
    ax.axvline(x=1 / 1000)
    plt.savefig("alpha_{}.pdf".format(alpha), bbox_inches="tight")
jacks = pd.read_csv("jks.csv", header=None)
jacks.columns = ["k = {}".format(i) for i in range(1,6)] + ["alpha", "T"]
dfs = []
for i in range(1, 6):
    cut = jacks[["k = {}".format(i), "alpha", "T"]]
    cut.columns = ["N_hat", "alpha", "T"]
    cut["estimator"] = "k = {}".format(i)
    dfs.append(cut)
dfs.append(estis.iloc[:,1:])
out = pd.concat(dfs)
out["1 / N_hat"] = 1 / out["N_hat"]
for alpha in out["alpha"].unique():
    fig, ax = plt.subplots(figsize=(18,12))
    sns.boxplot(x="1 / N_hat", y="T", hue="estimator", orient="h",
                data=out[out["alpha"] == alpha], whis=[5, 95], meanline=True,
                ax=ax)
    ax.axvline(x=1 / 1000)
    plt.savefig("jk_alpha_{}.pdf".format(alpha), bbox_inches="tight")
estis["log alpha_hat"] = np.log(estis["alpha_hat"])
for alpha in estis["alpha"].unique():
    fig, ax = plt.subplots(figsize=(18,12))
    sns.boxplot(x="log alpha_hat", y="T", orient="h",
                data=estis[estis["alpha"] == alpha], whis=[5, 95], meanline=True,
                ax=ax)
    #ax.axvline(x=alpha)
    plt.savefig("alpha_hat_{}.pdf".format(alpha), bbox_inches="tight")
    plt.close()
sns.boxplot()