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

DIR = "./_900_output/data/eqp_diffn/{}"
FIGS = "./_900_output/figures/eqp_diffn/{}"


def nan_trials(df):
    cut = df[df["N_hat"].isnull()][["T", "trial"]].to_numpy()
    return set(tuple(i) for i in cut)

def remove_indices(df, multindex, blacklist):
    df = df.set_index(multindex)
    return df.loc[set(df.index) - blacklist]    

sns.set()
chaos = pd.read_csv(DIR.format("chaos.csv"), header=None)
chaos.columns = ["N_hat", "T", "trial"]
chaos["estimator"] = "Chao"
chaos.replace(np.inf, value=np.nan, inplace=True)
fails = nan_trials(chaos)
chaos_corr = pd.read_csv(DIR.format("chaos_corr.csv"), header=None)
chaos_corr.columns = ["N_hat", "T", "trial"]
chaos_corr["estimator"] = "Chao (corrected)"
fails.update(nan_trials(chaos_corr))
links = pd.read_csv(DIR.format("links.csv"), header=None)
links.columns = ["N_hat", "T", "trial"]
links["estimator"] = "Lincoln & Schnabel"
links.replace(np.inf, value=np.nan, inplace=True)
fails.update(nan_trials(links))
schnab = pd.read_csv(DIR.format("schnab.csv"), header=None)
schnab.columns = ["N_hat", "T", "trial"]
schnab["estimator"] = "Lincoln & Schnabel"
schnab.replace(np.inf, value=np.nan, inplace=True)
fails.update(nan_trials(schnab))
estis = pd.read_csv(DIR.format("estis.csv"), header=None)
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
fig.savefig(FIGS.format("estimates_boxplots.pdf"), bbox_inches="tight")
plt.close(fig)
jacks = pd.read_csv(DIR.format("jks.csv"), header=None)
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
fig.savefig(FIGS.format("jks_boxplots.pdf"), bbox_inches="tight")
plt.close(fig)
estis["log alpha_hat"] = np.log(estis["alpha_hat"])
fig, ax = plt.subplots(figsize=(18,12))
sns.boxplot(x="log alpha_hat", y="T", orient="h",
            data=estis, whis=[5, 95], meanline=True,
            ax=ax)
ax.axvline(x=np.log(10.0))
fig.savefig(FIGS.format("alpha_hat.pdf"), bbox_inches="tight")
plt.close(fig)
nu_trace = data = pd.read_csv(DIR.format("Nu_trace.csv"), header=None)
nu_trace.columns = ([str(i) for i in list(range(0, 5100, 100))]
                    + ["alpha", "Nu", "O", "trial"])
nu_trace["T"] = [2, 4, 10] * 100
nu_trace["true_Nu"] = 1000 - nu_trace["O"]
for T in [2, 4, 10]:
    fig, ax = plt.subplots(figsize=(18,12))
    cut = nu_trace[nu_trace["T"] == T]
    for row in range(cut.shape[0]):
        ax.plot(range(0, 5100, 100), cut.iloc[row, :51])
        ax.axvline(x=cut["true_Nu"].iloc[i], color="gold", alpha=0.1)
        ax.set_ylabel("Loglikelihood")
        ax.set_xlabel(r"$N_U$")
        ax.set_title(r"Loglikelihood for different values of $N_U$ at $\alpha = \hat{\alpha}_{ML}$ for each trial. Vertical line shows true $N_U$")
    fig.savefig(FIGS.format("Nu_trace_{}.pdf".format(T)), bbox_inches="tight")
    plt.close(fig)
alpha_trace = pd.read_csv(DIR.format("alpha_trace.csv"), header=None)
alpha_trace.columns = ([str(i) for i in list(np.arange(0.1, 201.1, 1))] 
                       + ["alpha", "Nu", "O", "trial"])
alpha_trace["T"] = [2, 4, 10] * 100
for T in [2, 4, 10]:
    fig, ax = plt.subplots(figsize=(18,12))
    cut = alpha_trace[alpha_trace["T"] == T]
    for row in range(cut.shape[0]):
        ax.plot(np.arange(1.1, 201.1, 1), cut.iloc[row, 1:201])
        ax.set_ylabel("Loglikelihood")
        ax.set_xlabel(r"$\alpha$")
        ax.set_title(r"Loglikelihood for different values of $\alpha$ at $N_U = \hat{N}_{U(ML)}$ for each trial.")
    fig.savefig(FIGS.format("alpha_trace_{}.pdf".format(T)), bbox_inches="tight")
    plt.close(fig)
