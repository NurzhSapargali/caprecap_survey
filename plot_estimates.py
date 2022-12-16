#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 15:43:55 2022

@author: nurzhan
"""
from glob import glob

from matplotlib import rc

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

ESTIMATES = "./_900_output/data/"
N = 1000
ALPHA_RANGE = np.arange(0.1, 51.1)
NU_RANGE = np.arange(0, 5200, 200)
FIGSIZE = (18, 12)
MODELS = ["eqp", "diffp"]
LATEX_MAP = {"Pseudolikelihood": r"$\widehat{N}$",
             "Schnabel": r"$\widehat{N}_\mathrm{Schnab.}$",
             "Chao": r"$\widehat{N}_\mathrm{Chao}$",
             "Zelterman": r"$\widehat{N}_\mathrm{Zelt.}$",
             "Conway-Maxwell-Poisson": r"$\widehat{N}^*_\mathrm{LCMP}$",
             "Huggins": r"$\widehat{N}_\mathrm{Hugg.}$",
             "Turing Geometric": r"$\widehat{N}_\mathrm{TG}$",
             "Turing": r"$\widehat{N}_\mathrm{TP}$",
             "Jackknife k = 1": r"$\widehat{N}_\mathrm{JK,1}$",
             "Jackknife k = 2": r"$\widehat{N}_\mathrm{JK,2}$",
             "Jackknife k = 3": r"$\widehat{N}_\mathrm{JK,3}$",
             "Jackknife k = 4": r"$\widehat{N}_\mathrm{JK,4}$",
             "Jackknife k = 5": r"$\widehat{N}_\mathrm{JK,5}$"}


def plot_estimates(figsize, data, meanprops, true_pop, true_alpha):
    fig, ax = plt.subplots(figsize=figsize)
    sns.boxplot(x="1 / N_hat",
                y="T",
                orient="h",
                data=data,
                hue="label",
                showmeans=True,
                ax=ax,
                meanprops=meanprops,
                whis=[5, 95])
    ax.invert_yaxis()
    ax.axvline(x = 1.0 / true_pop, label=r"$N^{-1} = " + str(1.0 / true_pop) + "$",
               color="black", linestyle="dashed")
    ax.set_xlabel(r"$\widehat{N}^{-1}$")
    if model == "eqp":
        ax.set_title(r"$\alpha \to \infty$")
    else:
        ax.set_title(r"$\alpha = {}$".format(true_alpha))
    ax.legend()

def plot_alphas(figsize, data, meanprops, true_alpha):
    fig, ax = plt.subplots(figsize=figsize)
    meanprops["markersize"] = 10
    sns.boxplot(x="log alpha_hat",
                y="T",
                orient="h",
                data=data,
                showmeans=True,
                meanprops=meanprops,
                whis=[5, 95],
                ax=ax)
    ax.set_xlabel(r"$\log(\widehat{\alpha})$")
    ax.invert_yaxis()
    ax.axvline(x = np.log(true_alpha), color="black", linestyle="dashed")
    if model == "eqp":
        ax.set_title(r"$\alpha \to \infty$")
    else:
        ax.set_title(r"$\alpha = {}$".format(true_alpha))

meanprops = dict(marker="X", markerfacecolor="#4c4c4c", markeredgecolor="#4c4c4c")
sns.set(context="paper", style="ticks", font_scale=2.25)
rc('font',**{'family':'STIXGeneral'})
rc('mathtext',**{'fontset':'stix'})
summary_tables = {}
for model in MODELS:
    folder = ESTIMATES + model + "/"
    bench_colnames = ["N_hat", "N_o", "T"]
    if model != "eqp":
        bench_colnames.append("alpha")
    bench_colnames.append("estimator")
    bench = pd.concat([pd.read_csv(i, header=None, names=bench_colnames)
                       for i in glob(folder + "benchmarks*")])
    estis_colnames = ["alpha_hat", "N_u", "N_o", "T"]
    if model != "eqp":
        estis_colnames.append("alpha")
    estis = pd.concat([pd.read_csv(i, header=None, names=estis_colnames)
                       for i in glob(folder + "estimates*")])
    estis["N_hat"] = estis["N_u"] + estis["N_o"]
    estis["estimator"] = "Pseudolikelihood"
    df = pd.concat([estis, bench])
    df["N_hat"] = pd.to_numeric(df["N_hat"], errors="coerce")
    df = df.astype({"T": "int64"})
    df["label"] = df["estimator"].map(LATEX_MAP)
    if model == "eqp":
        df["alpha"] = np.inf
    df = df.iloc[::-1, :]
    for alpha in df["alpha"].unique():
        cut = df[df["alpha"] == alpha]
        print(alpha, cut[(cut["N_hat"].isna()) & (cut["T"] == 20)].groupby("estimator").count())
        cut["1 / N_hat"] = 1.0 / cut["N_hat"]
        cut["log alpha_hat"] = np.log(cut["alpha_hat"])
        plot_alphas(FIGSIZE, cut, meanprops, alpha)
        plt.savefig(folder.replace("data", "figures") + "alpha_boxplots_{}.pdf".format(alpha),
                    bbox_inches="tight")
        plt.close()
        cut["sq.error.N_hat"] = (N - cut["N_hat"])**2
        cut["sq.error.alpha_hat"] = (alpha - cut["alpha_hat"])**2
        mean_table = cut.groupby(["estimator", "T"]).mean()[["N_hat", "alpha_hat", "sq.error.N_hat", "sq.error.alpha_hat"]]
        mean_table["rmse.N_hat"] = mean_table["sq.error.N_hat"] ** 0.5
        mean_table["rel.bias.N_hat"] = (mean_table["N_hat"] - N) / N
        mean_table["rmse.alpha_hat"] = mean_table["sq.error.alpha_hat"] ** 0.5
        mean_table["rel.bias.alpha_hat"] = (mean_table["alpha_hat"] - alpha) / alpha
        summary_tables[alpha] = mean_table
        jks = cut[(cut["estimator"].str.contains("Jackknife"))|(cut["estimator"] == "Pseudolikelihood")]
        cut = cut[~cut["estimator"].str.contains("Jackknife")]
        meanprops["markersize"] = 6
        plot_estimates(FIGSIZE, cut, meanprops, N, alpha)
        plt.savefig(folder.replace("data", "figures") + "estimate_boxplots_{}.pdf".format(alpha),
                    bbox_inches="tight")
        plt.close()
        plot_estimates(FIGSIZE, jks, meanprops, N, alpha)
        plt.savefig(folder.replace("data", "figures") + "jackknife_boxplots_{}.pdf".format(alpha),
                    bbox_inches="tight")
        plt.close()
for alpha in summary_tables:
    cut = summary_tables[alpha].reset_index()
    for T in cut["T"].unique():
        base = cut.loc[(cut["estimator"] == "Pseudolikelihood")&(cut["T"] == T), "rmse.N_hat"].iloc[0]
        cut.loc[cut["T"] == T, "log.rel.rmse"] = np.log(cut.loc[cut["T"] == T, "rmse.N_hat"] / base)
    fig, ax = plt.subplots(figsize=FIGSIZE)
    sns.lineplot(y="log.rel.rmse",
                 x="T",
                 hue="estimator",
                 data=cut[cut["estimator"] != "Pseudolikelihood"],
                 marker="o",
                 style="estimator")
    ax.set_xticks(cut["T"].unique())
    ax.axhline(y=0, color="black")
    ax.set_ylim(-2.5, 2.5)
    plt.savefig(folder.replace("data", "figures") + "estimate_rmse_{}.pdf".format(alpha),
                    bbox_inches="tight")
    plt.close()
        # if model == "eqp":
        #     nu_name = "Nu_trace.csv"
        #     alpha_name = "alpha_trace.csv"
        # else:
        #     nu_name = "Nu_trace_{}.csv".format(alpha)
        #     alpha_name = "alpha_trace_{}.csv".format(alpha)           
        # nu_trace = pd.read_csv(folder + nu_name,
        #                        header=None)
        # nu_trace.columns = list(NU_RANGE) + ["alpha_hat", "N_u", "N_o", "T", "alpha"]
        # nu_trace["true N_u"] = N - nu_trace["N_o"]
        # for t in nu_trace["T"].unique():
        #     draws = nu_trace[nu_trace["T"] == t]
        #     fig, ax = plt.subplots(figsize=FIGSIZE)
        #     for row in range(draws.shape[0]):
        #         ax.plot(NU_RANGE, draws.iloc[row, :len(NU_RANGE)])
        #         ax.axvline(x=draws["true N_u"].iloc[row], color="gold", alpha=0.1)
        #         ax.set_ylabel("Loglikelihood")
        #         ax.set_xlabel(r"$N_u$")
        #         ax.set_title(r"$\log L(\widehat{\alpha}, N_u),$" + r" $\alpha = {}, T = {}.$".format(alpha, t))
        #     fig.savefig(folder.replace("data", "figures") + "Nu_trace_{}_{}.pdf".format(t, alpha),
        #                 bbox_inches="tight")
        #     plt.close(fig)
        # alpha_trace = pd.read_csv(folder + alpha_name,
        #                           header=None)
        # alpha_trace.columns = list(ALPHA_RANGE) + ["alpha_hat", "N_u", "N_o", "T", "alpha"]
        # for t in alpha_trace["T"].unique():
        #     draws = alpha_trace[nu_trace["T"] == t]
        #     fig, ax = plt.subplots(figsize=FIGSIZE)
        #     for row in range(draws.shape[0]):
        #         ax.plot(ALPHA_RANGE, draws.iloc[row, :len(ALPHA_RANGE)])
        #         ax.axvline(x=alpha, color="gold", alpha=0.1)
        #         ax.set_ylabel("Loglikelihood")
        #         ax.set_xlabel(r"$\alpha$")
        #         ax.set_title(r"$\log L(\alpha, \widehat{N_u}),$" + r" $\alpha = {}, T = {}.$".format(alpha, t))
        #     fig.savefig(folder.replace("data", "figures") + "alpha_trace_{}_{}.pdf".format(t, alpha),
        #                 bbox_inches="tight")
        #     plt.close(fig)
