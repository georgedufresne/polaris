#!/usr/bin/env python3
"""analysis/ks_vs_age_32k_vs_4m.png: best 2D KS at each age with observed-sized (32k) and 4M-star model sCMDs,
for production (zero-point map on) and for Ying's E box without the map, each with its own bootstrap 3 sigma
line. Run from m15_scmd4m/ after both 4M runs and their bootstraps."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt     # noqa: E402
import pandas as pd                 # noqa: E402

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
H = f"{B}/m15_scmd4m"
lim = {int(k): v for k, v in json.load(open("/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json")).items()}
ORANGE, BLUE, INK, MUTED, GRID, BG = "#eb6834", "#2a78d6", "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": BG, "axes.facecolor": BG})


def curve(d):
    g = pd.read_csv(f"{d}/combined_2d_ks_results.csv",
                    usecols=["simulation_id", "age_gyr", "distance_modulus", "extinction_ebv", "ks_statistic"])
    g = g[g.age_gyr < g.simulation_id.map(lim) - 1e-6]
    boot = pd.read_csv(f"{d}/bootstrap_2d_ks_results.csv").ks_statistic
    return g.groupby("age_gyr").ks_statistic.min(), boot.mean() + 3 * boot.std()


fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
for ax, (title, small, big) in zip(axes, (("Production (zero-point map on)", f"{B}/m15_final_refined", f"{H}/prod"),
                                          ("Ying's E box, no zero-point map", f"{B}/m15_yingbox_nozp_refined", f"{H}/nozp"))):
    for d, col, name in ((small, ORANGE, "32k-star sCMDs"), (big, BLUE, "4M-star sCMDs")):
        c, cut = curve(d)
        c = c[c.index >= 10.5]
        ax.plot(c.index, c.values, color=col, lw=2, label=name)
        ax.axhline(cut, color=col, lw=1, ls="--")
        ax.annotate(f"3σ, {name.split()[0]}", (10.55, cut), xytext=(0, 3), textcoords="offset points",
                    color=MUTED, fontsize=8)
        i = c.idxmin()
        ax.plot([i], [c[i]], "o", ms=8, color=col, mec=BG, mew=2)
        ax.annotate(f"{c[i]:.4f} at {i:.1f} Gyr", (i, c[i]), xytext=(0, -15), textcoords="offset points",
                    color=INK, fontsize=8.5, ha="center")
    ax.set_title(title, color=INK, fontsize=10.5, loc="left")
    ax.set_xlabel("age (Gyr)")
    ax.grid(axis="y", color=GRID, lw=0.8)
axes[0].set_ylabel("best 2D KS at this age")
axes[0].set_ylim(0, 0.045)
axes[1].legend(frameon=False, loc="upper right")
fig.suptitle("M15: best fit by age, model sCMDs the size of the data (32k) vs 4 million stars", color=INK,
             fontsize=11, x=0.01, ha="left")
fig.tight_layout()
os.makedirs(f"{H}/analysis", exist_ok=True)
fig.savefig(f"{H}/analysis/ks_vs_age_32k_vs_4m.png", dpi=150)
print(f"{H}/analysis/ks_vs_age_32k_vs_4m.png")
