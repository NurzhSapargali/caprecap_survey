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

DIR = "./_900_output/data/diffp_eqn/{}"
FIGS = "./_900_output/figures/diffp_eqn/{}"
N = 3000
T = [5, 10, 15, 20]
ALPHAS = [0.1, 3.0, 10.0]
ALPHA_RANGE = np.arange(0.1, 29.1)
NU_RANGE = np.arange(0, 5200, 200)
K = range(1, 6)
REPLICATIONS = 100

def nan_trials(df):
    cut = df[df["N_hat"].isnull()][["T", "trial"]].to_numpy()
    return set(tuple(i) for i in cut)

def remove_indices(df, multindex, blacklist):
    df = df.set_index(multindex)
    return df.loc[set(df.index) - blacklist]    


sns.set()
for alpha in ALPHAS:
    chaos = pd.read_csv(DIR.format("chaos_{}.csv".format(alpha)),
                        header=None)
    chaos.columns = ["N_hat", "T", "trial"]
    chaos["estimator"] = "Chao"
    chaos.replace(np.inf, value=np.nan, inplace=True)
    fails = nan_trials(chaos)
    chaos_corr = pd.read_csv(DIR.format("chaos_corr_{}.csv".format(alpha)),
                             header=None)
    chaos_corr.columns = ["N_hat", "No", "T", "trial"]
    chaos_corr["estimator"] = "Chao (corrected)"
    fails.update(nan_trials(chaos_corr))
    links = pd.read_csv(DIR.format("links_{}.csv".format(alpha)),
                        header=None)
    links.columns = ["N_hat", "No", "No", "T", "trial"]
    links["estimator"] = "Lincoln & Schnabel"
    links.replace(np.inf, value=np.nan, inplace=True)
    fails.update(nan_trials(links))
    schnab = pd.read_csv(DIR.format("schnab_{}.csv".format(alpha)),
                         header=None)
    schnab.columns = ["N_hat", "No", "T", "trial"]
    schnab["estimator"] = "Lincoln & Schnabel"
    schnab.replace(np.inf, value=np.nan, inplace=True)
    fails.update(nan_trials(schnab))
    estis = pd.read_csv(DIR.format("estis_{}.csv".format(alpha)),
                        header=None)
    estis.columns = ["alpha_hat", "Nu_hat", "No", "T", "trial"]
    estis["N_hat"] = estis["Nu_hat"] + estis["No"]
    estis["estimator"] = "Pseudo-likelihood"
    dfs = [links, schnab, chaos, chaos_corr, estis]
    out = pd.concat([remove_indices(i, ["T", "trial"], fails) for i in dfs])
    out["1 / N_hat"] = 1.0 / out["N_hat"]
    out = out.reset_index()
    fig, ax = plt.subplots(figsize=(18,12))
    sns.boxplot(x="1 / N_hat", y="T", hue="estimator", orient="h",
                data=out, whis=[5, 95], meanline=True,
                ax=ax)
    ax.axvline(x=1.0 / N)
    fig.savefig(FIGS.format("estimates_boxplots_{}.pdf".format(alpha)),
                bbox_inches="tight")
    plt.close(fig)
    jacks = pd.read_csv(DIR.format("jks_{}.csv".format(alpha)),
                        header=None)
    jacks.columns = ["k = {}".format(i) for i in K] + ["No", "T", "trial"]
    dfs = []
    for k in K:
        cut = jacks[["k = {}".format(k), "T", "trial"]]
        cut.columns = ["N_hat", "T", "trial"]
        cut["estimator"] = "k = {}".format(k)
        dfs.append(cut)
    dfs.append(estis)
    out = pd.concat([remove_indices(i, ["T", "trial"], fails) for i in dfs])
    out["1 / N_hat"] = 1.0 / out["N_hat"]
    out = out.reset_index()
    fig, ax = plt.subplots(figsize=(18,12))
    sns.boxplot(x="1 / N_hat", y="T", hue="estimator", orient="h",
                data=out, whis=[5, 95], meanline=True,
                ax=ax)
    ax.axvline(x=1.0 / N)
    fig.savefig(FIGS.format("jks_boxplots_{}.pdf".format(alpha)),
                bbox_inches="tight")
    plt.close(fig)
    estis["log alpha_hat"] = np.log(estis["alpha_hat"])
    fig, ax = plt.subplots(figsize=(18,12))
    sns.boxplot(x="log alpha_hat", y="T", orient="h",
                data=estis, whis=[5, 95], meanline=True,
                ax=ax)
    ax.axvline(x=np.log(alpha))
    fig.savefig(FIGS.format("alpha_hat_{}.pdf".format(alpha)),
                bbox_inches="tight")
    plt.close(fig)
    nu_trace = pd.read_csv(DIR.format("Nu_trace_{}.csv".format(alpha)),
                           header=None)
    nu_trace.columns = ([str(i) for i in list(NU_RANGE)]
                        + ["alpha", "Nu", "No", "trial"])
    nu_trace["true_Nu"] = N - nu_trace["No"]
    for t in T:
        fig, ax = plt.subplots(figsize=(18,12))
        cut = nu_trace[nu_trace["T"] == t]
        for row in range(cut.shape[0]):
            ax.plot(NU_RANGE[1:], cut.iloc[row, 1:len(NU_RANGE)]) # Discard first values to take logs
            ax.axvline(x=cut["true_Nu"].iloc[t], color="gold", alpha=0.1)
            ax.set_ylabel("Loglikelihood")
            ax.set_xlabel(r"$N_U$")
            ax.set_title(r"Loglikelihood for different values of $N_U$ at $\hat{\alpha} = \alpha$ for each trial. Vertical line shows true $N_U$," + " T = {}.".format(t))
        fig.savefig(FIGS.format("Nu_trace_{}_{}.pdf".format(t, alpha)),
                    bbox_inches="tight")
        plt.close(fig)
    alpha_trace = pd.read_csv(DIR.format("alpha_trace_{}.csv".format(alpha)),
                              header=None)
    alpha_trace.columns = ([str(i) for i in list(ALPHA_RANGE)] 
                           + ["alpha", "Nu", "No", "trial"])
    for t in T:
        fig, ax = plt.subplots(figsize=(18,12))
        cut = alpha_trace[alpha_trace["T"] == t]
        ax.axvline(x=alpha, color="gold")
        for row in range(cut.shape[0]):
            ax.plot(ALPHA_RANGE, cut.iloc[row, 1:len(ALPHA_RANGE) + 1])
            ax.set_ylabel("Loglikelihood")
            ax.set_xlabel(r"$\alpha$")
            ax.set_title(r"Loglikelihood for different values of $\alpha$ at $\hat{N}_U = N_{U}$ for each trial. Vertical line shows true $\alpha$," + " T = {}.".format(t))
        fig.savefig(FIGS.format("alpha_trace_{}_{}.pdf".format(t, alpha)),
                    bbox_inches="tight")
        plt.close(fig)  
