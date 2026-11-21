#!/usr/bin/env python3
"""analysis/ks_vs_age_our_vs_ying_isochrones.png: best 2D KS at each age for M15 with Ying's E box and no zero-point
map, our DSEP4 isochrones (../m15_yingbox_nozp_refined) vs Ying's published isochrones (../m15_yingiso_nozp_refined),
each with its own bootstrap 3 sigma line. Ying's 2D KS age (13.39) is marked."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt     # noqa: E402
import pandas as pd                 # noqa: E402

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
H = f"{B}/m15_yingiso_nozp"
ORANGE, BLUE, INK, MUTED, GRID, BG = "#eb6834", "#2a78d6", "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": BG, "axes.facecolor": BG})


def curve(d, limf):
    lim = {int(k): v for k, v in json.load(open(limf)).items()}
    g = pd.read_csv(f"{d}/combined_2d_ks_results.csv", usecols=["simulation_id", "age_gyr", "ks_statistic"])
    g = g[g.age_gyr < g.simulation_id.map(lim) - 1e-6]
    boot = pd.read_csv(f"{d}/bootstrap_2d_ks_results.csv").ks_statistic
    return g.groupby("age_gyr").ks_statistic.min(), boot.mean() + 3 * boot.std()


fig, ax = plt.subplots(figsize=(8, 4.6))
cuts = []
for d, limf, col, name in ((f"{B}/m15_yingbox_nozp_refined", "/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json",
                            ORANGE, "our isochrones (DSEP4 rebuild)"),
                           (f"{B}/m15_yingiso_nozp_refined", "/dartfs-hpc/scratch/f007058/M15/iso_ying2025/age_limits_ying.json",
                            BLUE, "Ying's published isochrones")):
    c, cut = curve(d, limf)
    c = c[c.index >= 10.5]
    ax.plot(c.index, c.values, color=col, lw=2, label=name)
    ax.axhline(cut, color=col, lw=1, ls="--")
    cuts.append(cut)
    i = c.idxmin()
    ax.plot([i], [c[i]], "o", ms=8, color=col, mec=BG, mew=2)
    ax.annotate(f"{c[i]:.4f} at {i:.1f} Gyr", (i, c[i]), xytext=(0, 9 if col == ORANGE else -15), textcoords="offset points", color=INK,
                fontsize=8.5, ha="center")
if abs(cuts[0] - cuts[1]) < 0.0004:      # the two nulls coincide: one label
    ax.annotate(f"bootstrap 3σ line, both runs ({cuts[0]:.4f} / {cuts[1]:.4f})", (10.55, max(cuts)), xytext=(0, 3),
                textcoords="offset points", color=MUTED, fontsize=8)
else:
    for cut, nm, dy in zip(cuts, ("ours", "Ying's"), (3, -11)):
        ax.annotate(f"3σ ({nm})", (10.55, cut), xytext=(0, dy), textcoords="offset points", color=MUTED, fontsize=8)
ax.axvline(13.39, color=MUTED, lw=1, ls=":")
ax.annotate("Ying 2025, 2D KS: 13.39 Gyr", (13.39, 0.0435), xytext=(-4, 0), textcoords="offset points", color=MUTED,
            fontsize=8, ha="right")
ax.set_xlabel("age (Gyr)")
ax.set_ylabel("best 2D KS at this age")
ax.set_ylim(0, 0.045)
ax.grid(axis="y", color=GRID, lw=0.8)
ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 0.93))
ax.set_title("M15, Ying's E box, no zero-point map: our isochrones vs Ying's (same sCMDs and KS)", color=INK,
             fontsize=10.5, loc="left")
fig.tight_layout()
os.makedirs(f"{H}/analysis", exist_ok=True)
fig.savefig(f"{H}/analysis/ks_vs_age_our_vs_ying_isochrones.png", dpi=150)
print(f"{H}/analysis/ks_vs_age_our_vs_ying_isochrones.png")
