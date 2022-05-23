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


def nan_trials(df):
    cut = df[df["N_hat"].isnull()][["T", "trial"]].to_numpy()
    return set(tuple(i) for i in cut)

def remove_indices(df, multindex, blacklist):
    df = df.set_index(multindex)
    return df.loc[set(df.index) - blacklist]    

chaos = pd.read_csv("chaos.csv", header=None)
chaos.columns = ["N_hat", "T", "trial"]
chaos["estimator"] = "Chao"
chaos.replace(np.inf, value=np.nan, inplace=True)
fails = nan_trials(chaos)
chaos_corr = pd.read_csv("chaos_corr.csv", header=None)
chaos_corr.columns = ["N_hat", "T", "trial"]
chaos_corr["estimator"] = "Chao (corrected)"
fails.update(nan_trials(chaos_corr))
links = pd.read_csv("links.csv", header=None)
links.columns = ["N_hat", "T", "trial"]
links["estimator"] = "Lincoln & Schnabel"
links.replace(np.inf, value=np.nan, inplace=True)
fails.update(nan_trials(links))
schnab = pd.read_csv("schnab.csv", header=None)
schnab.columns = ["N_hat", "T", "trial"]
schnab["estimator"] = "Lincoln & Schnabel"
schnab.replace(np.inf, value=np.nan, inplace=True)
fails.update(nan_trials(schnab))
estis = pd.read_csv("estis.csv", header=None)
estis.columns = ["alpha_hat", "N_hat", "T", "trial"]
estis["estimator"] = "Pseudo-likelihood"
dfs = [links, schnab, chaos, chaos_corr, estis]
out = pd.concat([remove_indices(i, ["T", "trial"], fails) for i in dfs])
out["1 / N_hat"] = 1 / out["N_hat"]
out = out.reset_index()
fig, ax = plt.subplots(figsize=(18,12))
sns.boxplot(x="1 / N_hat", y="T", hue="estimator", orient="h",
                data=out, whis=[5, 95], meanline=True,
                ax=ax)
ax.axvline(x=1 / 1000)
fig.savefig("estimates_boxplots.pdf", bbox_inches="tight")
jacks = pd.read_csv("jks.csv", header=None)
jacks.columns = ["k = {}".format(i) for i in range(1,6)] + ["T", "trial"]
dfs = []
for i in range(1, 6):
    cut = jacks[["k = {}".format(i), "T", "trial"]]
    cut.columns = ["N_hat", "T", "trial"]
    cut["estimator"] = "k = {}".format(i)
    dfs.append(cut)
dfs.append(estis)
out = pd.concat([remove_indices(i, ["T", "trial"], fails) for i in dfs])
out["1 / N_hat"] = 1 / out["N_hat"]
out = out.reset_index()
fig, ax = plt.subplots(figsize=(18,12))
sns.boxplot(x="1 / N_hat", y="T", hue="estimator", orient="h",
                data=out, whis=[5, 95], meanline=True,
                ax=ax)
ax.axvline(x=1 / 1000)
fig.savefig("jks_boxplots.pdf", bbox_inches="tight")
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