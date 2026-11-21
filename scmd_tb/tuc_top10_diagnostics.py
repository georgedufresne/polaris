#!/usr/bin/env python3
"""Take the 10 best 47 Tuc runs apart: replay check, luminosity function, colour spread, where the 2D KS comes from."""
import json
import sys
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp

G = "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02"
OUT = "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_top10_diagnostics"
sys.path.insert(0, f"{G}/code")                       # the grid's frozen code: exact replays
from isochronetoolbox.scmd import Setup, context, regenerate   # noqa: E402
from isochronetoolbox.scmd.pipeline import _window              # noqa: E402
from isochronetoolbox.scmd.stats import ks_2d                   # noqa: E402

S = Setup(**{k: v for k, v in json.load(open(f"{G}/setup.json")).items() if k in Setup.__dataclass_fields__})
_ctx = None


def ctx():
    global _ctx
    if _ctx is None:
        _ctx = context(S)
    return _ctx


def one(r):
    c = ctx()
    cat = regenerate(c, int(r["simulation_id"]), r["age_gyr"], r["distance_modulus"], r["extinction_ebv"])
    sim = c.points(cat)
    sim = sim[_window(sim, S.ks_window)]
    return r, sim


def ridge_colour(ridge, mag, lo, hi):
    m = (ridge.mag > lo) & (ridge.mag < hi)
    o = np.argsort(ridge.mag[m])
    return np.interp(mag, ridge.mag[m][o], ridge.color[m][o])


def main():
    import os
    os.makedirs(OUT, exist_ok=True)
    comb = pd.read_csv(f"{G}/combined_2d_ks_results.csv")
    best = comb.loc[comb.groupby("simulation_id").ks_statistic.idxmin()].nsmallest(10, "ks_statistic")
    with Pool(5) as p:
        res = p.map(one, best.to_dict("records"))
    from types import SimpleNamespace
    res = [(SimpleNamespace(**r), sim) for r, sim in res]
    c = ctx()
    obs = c.obs3
    ridge = c.rcut.ridge
    ms_o, rgb_o = obs[obs[:, 1] > 17.5], obs[obs[:, 1] < 16.6]
    rows, sims = [], []
    for r, sim in res:
        ks, d = ks_2d(obs, sim, c.obs3_cdf)
        ks_mag = ks_2samp(obs[:, 1], sim[:, 1]).statistic                       # luminosity function
        ms_s, rgb_s = sim[sim[:, 1] > 17.5], sim[sim[:, 1] < 16.6]
        dc_ms_o = ms_o[:, 0] - ridge_colour(ridge, ms_o[:, 1], 17.4, 19.1)
        dc_ms_s = ms_s[:, 0] - ridge_colour(ridge, ms_s[:, 1], 17.4, 19.1)
        dc_rg_o = rgb_o[:, 0] - ridge_colour(ridge, rgb_o[:, 1], 14.9, 16.7)
        dc_rg_s = rgb_s[:, 0] - ridge_colour(ridge, rgb_s[:, 1], 14.9, 16.7)
        mad = lambda x: 1.4826 * np.median(np.abs(x - np.median(x)))
        rows.append({"run": int(r.simulation_id), "age": r.age_gyr, "mu": r.distance_modulus, "E(V-I)": r.extinction_ebv,
                     "KS table": r.ks_statistic, "KS replay": ks,
                     "KS F606W only (LF)": ks_mag,
                     "KS MS colour": ks_2samp(dc_ms_o, dc_ms_s).statistic,
                     "KS RGB colour": ks_2samp(dc_rg_o, dc_rg_s).statistic,
                     "MS colour offset": np.median(dc_ms_s) - np.median(dc_ms_o),
                     "MS width ratio": mad(dc_ms_s) / mad(dc_ms_o),
                     "RGB colour offset": np.median(dc_rg_s) - np.median(dc_rg_o),
                     "N(V<17) model/obs": np.mean(sim[:, 1] < 17.0) / np.mean(obs[:, 1] < 17.0)})
        sims.append((r, sim, d))
    tab = pd.DataFrame(rows)
    tab.to_csv(f"{OUT}/top10_diagnostics.csv", index=False)
    pd.set_option("display.width", 250)
    print(tab.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    # figure: LF ratio, MS colour residuals, where the KS differences sit (best 3)
    fig, ax = plt.subplots(2, 3, figsize=(17, 10))
    edges = np.arange(15.0, 19.05, 0.1)
    ho, _ = np.histogram(obs[:, 1], edges)
    cen = 0.5 * (edges[1:] + edges[:-1])
    for r, sim, d in sims[:5]:
        hs, _ = np.histogram(sim[:, 1], edges)
        ax[0, 0].plot(cen, hs / np.maximum(ho, 1), lw=1.2, label=f"run {int(r.simulation_id)} {r.age_gyr} Gyr")
    ax[0, 0].axhline(1, color="0.3", lw=0.8)
    ax[0, 0].axvline(17.02, color="0.6", ls=":", lw=1)
    ax[0, 0].set_xlabel("F606W")
    ax[0, 0].set_ylabel("model / observed stars per 0.1 mag")
    ax[0, 0].set_title("luminosity function: model over data (5 best runs)\n(dotted: SGB)", fontsize=10)
    ax[0, 0].set_ylim(0.6, 1.6)
    ax[0, 0].legend(fontsize=7)
    b = np.linspace(-0.08, 0.08, 81)
    ax[0, 1].hist(obs[obs[:, 1] > 17.5, 0] - ridge_colour(ridge, obs[obs[:, 1] > 17.5, 1], 17.4, 19.1), b, histtype="step",
                  color="k", lw=1.8, density=True, label="observed")
    for r, sim, d in sims[:3]:
        s = sim[sim[:, 1] > 17.5]
        ax[0, 1].hist(s[:, 0] - ridge_colour(ridge, s[:, 1], 17.4, 19.1), b, histtype="step", lw=1.2, density=True,
                      label=f"run {int(r.simulation_id)}")
    ax[0, 1].set_yscale("log")
    ax[0, 1].set_xlabel("colour - observed ridge (main sequence, F606W > 17.5)")
    ax[0, 1].set_title("main-sequence colour spread", fontsize=10)
    ax[0, 1].legend(fontsize=8)
    ax[0, 2].axis("off")
    for a, (r, sim, d) in zip(ax[1], sims[:3]):
        o = np.argsort(d)
        sc = a.scatter(obs[o, 0], obs[o, 1], c=d[o], s=0.4, cmap="viridis", vmin=0, vmax=np.percentile(d, 99), rasterized=True)
        top = d >= np.percentile(d, 95)
        a.set_xlim(0.35, 0.95)
        a.set_ylim(19.05, 14.95)
        a.set_title(f"run {int(r.simulation_id)}: |ECDF difference| at each observed star\n"
                    f"KS = 95th percentile = {np.percentile(d, 95):.4f}; top 5% at F606W "
                    f"{np.percentile(obs[top, 1], 10):.2f}-{np.percentile(obs[top, 1], 90):.2f}", fontsize=9)
        a.set_xlabel("F606W - F814W")
    ax[1, 0].set_ylabel("F606W")
    fig.colorbar(sc, ax=ax[1].tolist(), shrink=0.8, label="|ECDF_obs - ECDF_model|")
    fig.savefig(f"{OUT}/top10_diagnostics.png", dpi=110, bbox_inches="tight")
    print("figure:", f"{OUT}/top10_diagnostics.png")


if __name__ == "__main__":
    main()
