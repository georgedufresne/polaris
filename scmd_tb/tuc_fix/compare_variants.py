#!/usr/bin/env python3
"""How far each AST-model fix moves the 47 Tuc fit (top-10 v3 runs, full grids).

Variants, each adding to the previous one:
  A_v3              v3 as it was (rerun with the same code, so differences are the fixes)
  B_zp              zero points of the fitted catalogue (.cal.adj.zpt, deep frames):
                    zp_v 30.5742 -> 30.5911, zp_i 29.6777 -> 29.7079
  C_zp_mix          no hard short/long-exposure split in the AST neighbour search
  D_zp_mix_nearest  plain resampling of one neighbour's paired errors (no kernel)

Writes variant_best.csv, variant_runs.csv and ks_vs_age.png here.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

F = os.path.dirname(os.path.abspath(__file__))
V = ["A_v3", "B_zp", "C_zp_mix", "D_zp_mix_nearest"]
LAB = {"A_v3": "A  v3", "B_zp": "B  + zero points", "C_zp_mix": "C  + no exposure split",
       "D_zp_mix_nearest": "D  + plain resampling"}
COL = {"A_v3": "#2a78d6", "B_zp": "#eb6834", "C_zp_mix": "#1baf7a", "D_zp_mix_nearest": "#eda100"}
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
KS3 = 0.0105   # v3 bootstrap mean + 3 sigma (setup ks_3sigma)


def load(v):
    p = os.path.join(F, v, "combined_2d_ks_results.csv")
    return pd.read_csv(p)


def main():
    tabs = {v: load(v) for v in V}
    rows, per = [], []
    for v, g in tabs.items():
        b = g.loc[g.ks_statistic.idxmin()]
        rows.append(dict(variant=v, run=int(b.simulation_id), age=b.age_gyr, dm=b.distance_modulus,
                         ebv=b.extinction_ebv, ks=b.ks_statistic, n_rows=len(g),
                         n_below_3sigma=int((g.ks_statistic <= KS3).sum())))
        for r, gr in g.groupby("simulation_id"):
            k = gr.ks_statistic.idxmin()
            per.append(dict(variant=v, run=r, age=gr.loc[k, "age_gyr"], dm=gr.loc[k, "distance_modulus"],
                            ebv=gr.loc[k, "extinction_ebv"], ks=gr.loc[k, "ks_statistic"]))
    best, per = pd.DataFrame(rows), pd.DataFrame(per)
    best.to_csv(os.path.join(F, "variant_best.csv"), index=False)
    per.to_csv(os.path.join(F, "variant_runs.csv"), index=False)
    pd.set_option("display.width", 200)
    print("best cell per variant")
    print(best.to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    wide = per.pivot(index="run", columns="variant")
    print("\nper run: best KS, and the age / distance modulus / E(B-V) of that cell")
    for q in ("ks", "age", "dm", "ebv"):
        t = wide[q][V].copy()
        t.columns = [c.split("_")[0] for c in V]
        print(f"\n{q}"); print(t.to_string(float_format=lambda x: f"{x:.4f}" if q in ('ks', 'ebv') else f"{x:.2f}"))
    for v in V[1:]:
        d = wide["age"][v] - wide["age"]["A_v3"]
        dk = wide["ks"][v] - wide["ks"]["A_v3"]
        dd = wide["dm"][v] - wide["dm"]["A_v3"]
        print(f"\n{v} vs A over the 10 runs: best-age shift median {d.median():+.2f} Gyr (range {d.min():+.1f}..{d.max():+.1f}),"
              f" DM shift median {dd.median():+.3f}, best-KS change median {dk.median():+.5f}")

    # KS profile: lowest KS over (DM, E) at each age, for the 4 best v3 runs
    top = per[per.variant == "A_v3"].sort_values("ks").run.head(4).tolist()
    plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
                         "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2,
                         "ytick.color": INK2, "text.color": INK, "axes.grid": True, "grid.color": GRID,
                         "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
                         "font.size": 9})
    fig, axs = plt.subplots(1, len(top), figsize=(3.4 * len(top), 3.2), sharey=True)
    for ax, r in zip(axs, top):
        for v in V:
            g = tabs[v][tabs[v].simulation_id == r].groupby("age_gyr").ks_statistic.min()
            ax.plot(g.index, g.values, "-", color=COL[v], lw=1.5, label=LAB[v])
        ax.axhline(KS3, color=INK2, lw=0.8, ls=":")
        ax.set_title(f"run {r}", loc="left")
        ax.set_xlabel("age (Gyr)")
    axs[0].set_ylabel("lowest 2D KS over distance and reddening")
    axs[0].annotate("v3 3$\\sigma$ cut", (axs[0].get_xlim()[0], KS3), xytext=(3, 3), textcoords="offset points",
                    color=INK2, fontsize=7)
    fig.legend(*axs[0].get_legend_handles_labels(), loc="upper right", ncol=4, frameon=False, fontsize=8)
    fig.suptitle("47 Tuc: effect of each AST-model fix on the fit", x=0.01, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(os.path.join(F, "ks_vs_age.png"), dpi=160)
    print(os.path.join(F, "ks_vs_age.png"))


if __name__ == "__main__":
    main()
