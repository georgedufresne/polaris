#!/usr/bin/env python3
"""analysis/sigma_vs_age.png: at each age, the best fit's distance above its own bootstrap null,
(best KS - null mean) / null std, for three M15 runs with Ying's E box and no zero-point map:
our isochrones (../m15_yingbox_nozp_refined), Ying's isochrones (../m15_yingiso_nozp_refined), and Ying's isochrones +
his artificial-star error model + 99th-percentile KS (../m15_yingiso_err99_refined). Sigma units put the 95th- and
99th-percentile statistics on one scale; the 3 sigma line is shared. Palette: categorical slots 1-3, validated
(validate_palette.js: all checks pass; aqua below 3:1 contrast, so its minimum is labelled directly)."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt     # noqa: E402
import pandas as pd                 # noqa: E402

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
H = f"{B}/m15_yingiso_err99"
D4 = "/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json"
YL = "/dartfs-hpc/scratch/f007058/M15/iso_ying2025/age_limits_ying.json"
ORANGE, BLUE, AQUA, INK, MUTED, GRID, BG = "#eb6834", "#2a78d6", "#1baf7a", "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": BG, "axes.facecolor": BG})
RUNS = [(f"{B}/m15_yingbox_nozp_refined", D4, ORANGE, "our isochrones, our errors, KS p95"),
        (f"{B}/m15_yingiso_nozp_refined", YL, BLUE, "Ying's isochrones, our errors, KS p95"),
        (f"{B}/m15_yingiso_err99_refined", YL, AQUA, "Ying's isochrones, his errors, KS p99")]

fig, ax = plt.subplots(figsize=(8.4, 4.8))
for d, limf, col, name in RUNS:
    lim = {int(k): v for k, v in json.load(open(limf)).items()}
    g = pd.read_csv(f"{d}/combined_2d_ks_results.csv", usecols=["simulation_id", "age_gyr", "ks_statistic"])
    g = g[g.age_gyr < g.simulation_id.map(lim) - 1e-6]
    boot = pd.read_csv(f"{d}/bootstrap_2d_ks_results.csv").ks_statistic
    c = (g.groupby("age_gyr").ks_statistic.min() - boot.mean()) / boot.std()
    c = c[c.index >= 10.5]
    ax.plot(c.index, c.values, color=col, lw=2, label=name)
    i = c.idxmin()
    ax.plot([i], [c[i]], "o", ms=8, color=col, mec=BG, mew=2)
    ax.annotate(f"{c[i]:.1f}σ at {i:.1f} Gyr", (i, c[i]), xytext=(0, 9 if col == ORANGE else -15), textcoords="offset points", color=INK,
                fontsize=8.5, ha="center")
ax.axhline(3, color=MUTED, lw=1, ls="--")
ax.annotate("3σ above the null (strict cut)", (14.5, 3), xytext=(0, 3), textcoords="offset points", color=MUTED, fontsize=8)
ax.axvline(13.39, color=MUTED, lw=1, ls=":")
ax.annotate("Ying 2025, 2D KS: 13.39 Gyr", (13.39, 11.6), xytext=(-4, 0), textcoords="offset points", color=MUTED,
            fontsize=8, ha="right")
ax.set_xlabel("age (Gyr)")
ax.set_ylabel("best fit at this age: σ above its own null")
ax.set_ylim(0, 12)
ax.grid(axis="y", color=GRID, lw=0.8)
ax.legend(frameon=False, loc="lower right")
ax.set_title("M15, Ying's E box, no zero-point map: isochrones, error model and KS percentile", color=INK,
             fontsize=10.5, loc="left")
fig.tight_layout()
os.makedirs(f"{H}/analysis", exist_ok=True)
fig.savefig(f"{H}/analysis/sigma_vs_age.png", dpi=150)
print(f"{H}/analysis/sigma_vs_age.png")
