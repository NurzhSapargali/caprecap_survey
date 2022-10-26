#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 15:43:55 2022

@author: nurzhan
"""
from glob import glob

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

ESTIMATES = "./_900_output/data/"
N = 1000
ALPHA_RANGE = np.arange(0.1, 51.1)
NU_RANGE = np.arange(0, 5200, 200)
FIGSIZE = (18, 12)

def nan_trials(df):
    cut = df[df["N_hat"].isnull()][["T", "trial"]].to_numpy()
    return set(tuple(i) for i in cut)

def remove_indices(df, multindex, blacklist):
    df = df.set_index(multindex)
    return df.loc[set(df.index) - blacklist]    


sns.set_theme(context="talk", style="whitegrid")
for folder in [ESTIMATES + "diffp/"]:
    print(folder)
    bench_colnames = ['N_hat', 'N_o', 'T', 'alpha', 'estimator']
    bench = pd.concat([pd.read_csv(i, header=None, names=bench_colnames)
                       for i in glob(folder + "benchmarks*.csv")])
    estis_colnames = ["alpha_hat", "N_u", "N_o", "T", "alpha"]
    estis = pd.concat([pd.read_csv(i, header=None, names=estis_colnames)
                       for i in glob(folder + "estimates*.csv")])
    estis["N_hat"] = estis["N_u"] + estis["N_o"]
    estis["estimator"] = "Pseudolikelihood"
    df = pd.concat([estis, bench])
    for alpha in df["alpha"].unique():
        cut = df[(df["alpha"] == alpha)&(df["estimator"] != "Conway-Maxwell-Poisson")]
        cut["N_hat"] = pd.to_numeric(cut["N_hat"], errors='coerce')
        cut["1 / N_hat"] = 1.0 / cut["N_hat"]
        cut["log alpha_hat"] = np.log(cut["alpha_hat"])
        fig, ax = plt.subplots(figsize=FIGSIZE)
        sns.boxplot(x="log alpha_hat",
                    y="T",
                    orient="h",
                    data=cut,
                    showmeans=True,
                    ax=ax,
                    meanprops={"markerfacecolor": "black"})
        ax.invert_yaxis()
        ax.axvline(x = np.log(alpha))
        fig.savefig(folder.replace("data", "figures") + "alpha_boxplots_{}.pdf".format(alpha),
                    bbox_inches="tight")
        plt.close(fig)
        jks = cut[(cut["estimator"].str.contains("Jackknife"))|(cut["estimator"] == "Pseudolikelihood")]
        cut = cut[~cut["estimator"].str.contains("Jackknife")]
        fig, ax = plt.subplots(figsize=FIGSIZE)
        sns.boxplot(x="1 / N_hat",
                    y="T",
                    orient="h",
                    data=cut,
                    hue="estimator",
                    showmeans=True,
                    ax=ax,
                    meanprops={"markerfacecolor": "black"})
        ax.invert_yaxis()
        ax.axvline(x = 1.0 / N)
        ax.set_title("Alpha = {}".format(alpha))
        fig.savefig(folder.replace("data", "figures") + "estimate_boxplots_{}.pdf".format(alpha),
                    bbox_inches="tight")
        plt.close(fig)
        fig, ax = plt.subplots(figsize=FIGSIZE)
        sns.boxplot(x="1 / N_hat",
                    y="T",
                    orient="h",
                    data=jks,
                    hue="estimator",
                    showmeans=True,
                    ax=ax,
                    meanprops={"markerfacecolor": "black"})
        ax.invert_yaxis()
        ax.axvline(x = 1.0 / N)
        ax.set_title("Alpha = {}".format(alpha))
        fig.savefig(folder.replace("data", "figures") + "jackknife_boxplots_{}.pdf".format(alpha),
                    bbox_inches="tight")
        plt.close(fig)
        nu_trace = pd.read_csv(folder + "Nu_trace_{}.csv".format(alpha),
                               header=None)
        nu_trace.columns = list(NU_RANGE) + ["alpha_hat", "N_u", "N_o", "T", "alpha"]
        nu_trace["true N_u"] = N - nu_trace["N_o"]
        for t in nu_trace["T"].unique():
            draws = nu_trace[nu_trace["T"] == t]
            fig, ax = plt.subplots(figsize=FIGSIZE)
            for row in range(draws.shape[0]):
                ax.plot(NU_RANGE, draws.iloc[row, :len(NU_RANGE)])
                ax.axvline(x=draws["true N_u"].iloc[row], color="gold", alpha=0.1)
                ax.set_ylabel("Loglikelihood")
                ax.set_xlabel(r"$N_u$")
                ax.set_title(r"$\log L(\widehat{\alpha}, N_u),$" + r" $\alpha = {}, T = {}.$".format(alpha, t))
            fig.savefig(folder.replace("data", "figures") + "Nu_trace_{}_{}.pdf".format(t, alpha),
                        bbox_inches="tight")
            plt.close(fig)
        alpha_trace = pd.read_csv(folder + "alpha_trace_{}.csv".format(alpha),
                                  header=None)
        alpha_trace.columns = list(ALPHA_RANGE) + ["alpha_hat", "N_u", "N_o", "T", "alpha"]
        for t in alpha_trace["T"].unique():
            draws = alpha_trace[nu_trace["T"] == t]
            fig, ax = plt.subplots(figsize=FIGSIZE)
            for row in range(draws.shape[0]):
                ax.plot(ALPHA_RANGE, draws.iloc[row, :len(ALPHA_RANGE)])
                ax.axvline(x=alpha, color="gold", alpha=0.1)
                ax.set_ylabel("Loglikelihood")
                ax.set_xlabel(r"$\alpha$")
                ax.set_title(r"$\log L(\alpha, \widehat{N_u}),$" + r" $\alpha = {}, T = {}.$".format(alpha, t))
            fig.savefig(folder.replace("data", "figures") + "alpha_trace_{}_{}.pdf".format(t, alpha),
                        bbox_inches="tight")
            plt.close(fig)
