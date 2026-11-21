#!/usr/bin/env python3
"""Figures for a refined sCMD grid (<grid>_final_refined/, from finish_refine_one.sh).

    refined_plots.py VIEW_DIR AGE_LIMITS_JSON LABEL

1. analysis/ks_vs_age.png: best KS at each age (minimum over runs, mu, E), stage 2 + 3 grid alone against
   the KS refinement, with the bootstrap null (mean, 3 sigma line); below, how many runs have a fit
   inside 3 sigma at each age (the cells that carry weight in the strict posterior).
2. analysis/passing_cells.png: age against E and against mu for every (run, age) minimum inside 3 sigma,
   shaded by KS, over all (run, age) minima in grey.
"""
import json
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

V, LIMF, LABEL = sys.argv[1].rstrip("/"), sys.argv[2], sys.argv[3]
BLUE, ORANGE, INK, MUTED, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#e4e3df"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": "#fcfcfb", "axes.facecolor": "#fcfcfb"})

lim = {int(k): v for k, v in json.load(open(LIMF)).items()}
m = pd.read_csv(f"{V}/combined_2d_ks_results.csv")
m = m[m.age_gyr < m.simulation_id.map(lim) - 1e-6]
boot = pd.read_csv(f"{V}/bootstrap_2d_ks_results.csv").ks_statistic.values
mu, sd = boot.mean(), boot.std()
cut = mu + 3 * sd


def per_cell(rows):
    rows = rows.reset_index(drop=True)
    return rows.loc[rows.groupby(["simulation_id", "age_gyr"]).ks_statistic.idxmin()]


best = per_cell(m)
grid = per_cell(m[m.source == "grid"]) if "source" in m else None
ok = best[best.ks_statistic <= cut]

# 1. KS against age, and runs inside 3 sigma per age
fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.2, 6.4), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
a1.axhspan(mu - sd, mu + sd, color=GRID, lw=0)
a1.axhline(mu, color=MUTED, lw=1)
a1.axhline(cut, color=MUTED, lw=1.2, ls="--")
if grid is not None:
    gb = grid.groupby("age_gyr").ks_statistic.min()
    a1.plot(gb.index, gb.values, color=ORANGE, lw=2, label="grid (stage 2 + 3)")
rb = best.groupby("age_gyr").ks_statistic.min()
a1.plot(rb.index, rb.values, color=BLUE, lw=2, label="refined")
i = rb.idxmin()
a1.plot([i], [rb[i]], "o", ms=8, color=BLUE, mec="#fcfcfb", mew=2)
a1.annotate(f"best {rb[i]:.4f} at {i:.1f} Gyr", (i, rb[i]), xytext=(0, -16), textcoords="offset points",
            color=INK, fontsize=9, ha="center")
lo = min(rb.min(), mu - sd) * 0.9
a1.set_ylim(lo, max(3 * cut, rb[rb.index > rb.idxmin() - 2].max() * 1.05))
a1.set_ylabel("best 2D KS at this age")
a1.legend(frameon=False, loc="upper right", fontsize=9)
a1.grid(axis="y", color=GRID, lw=0.6)
a1.set_title(f"{LABEL}: best fit by age against the bootstrap null ({len(boot):,} realisations)", fontsize=10.5,
             loc="left", color=INK)
n = ok.groupby("age_gyr").simulation_id.nunique().reindex(rb.index, fill_value=0)
a2.bar(n.index, n.values, width=0.08, color=BLUE, lw=0)
a2.set_ylabel("runs inside 3σ")
a2.set_xlabel("age (Gyr)")
a2.grid(axis="y", color=GRID, lw=0.6)
a2.text(0.01, 0.92, f"{len(ok)} (run, age) cells from {ok.simulation_id.nunique()} runs inside 3σ",
        transform=a2.transAxes, va="top", color=INK, fontsize=9)
span = ok.age_gyr if len(ok) else rb.index.to_series()
a2.set_xlim(max(rb.index.min(), span.min() - 1.5), min(rb.index.max(), span.max() + 1.5) + 0.05)
x0 = a2.get_xlim()[0]
a1.text(x0, mu + sd, " null mean ± 1σ", va="bottom", ha="left", color=MUTED, fontsize=9)
a1.text(x0, cut, f" null 3σ line, KS {cut:.4f}", va="top", ha="left", color=MUTED, fontsize=9)
fig.tight_layout()
fig.savefig(f"{V}/analysis/ks_vs_age.png", dpi=150)
plt.close(fig)

# 2. where the passing cells sit in age-E and age-mu
fig, axes = plt.subplots(1, 2, figsize=(10, 4.4), sharex=True, gridspec_kw={"wspace": 0.32})
near = best[best.ks_statistic <= 2 * cut]
for ax, col, lab in [(axes[0], "extinction_ebv", "E(F606W-F814W) = 1.02 E(B-V)"),
                     (axes[1], "distance_modulus", "apparent F606W modulus")]:
    jit = np.random.default_rng(1).uniform(-0.03, 0.03, len(near))
    ax.scatter(near.age_gyr + jit, near[col], s=6, color="#c3c2b7", lw=0, label="(run, age) minima, KS < 2× the 3σ line")
    o = ok.sort_values("ks_statistic", ascending=False)
    jo = np.random.default_rng(2).uniform(-0.03, 0.03, len(o))
    sc = ax.scatter(o.age_gyr + jo, o[col], c=o.ks_statistic, cmap="Blues_r", vmin=o.ks_statistic.min(), vmax=cut * 1.15,
                    s=16, edgecolors="#fcfcfb", linewidths=0.4, label="inside 3σ")
    ax.set_xlabel("age (Gyr)")
    ax.set_ylabel(lab)
    ax.grid(color=GRID, lw=0.6)
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, frameon=False, fontsize=8.5, loc="lower center", ncol=2, markerscale=1.5, bbox_to_anchor=(0.45, -0.06))
cb = fig.colorbar(sc, ax=axes, shrink=0.85, pad=0.015)
cb.set_label("2D KS")
if len(ok):
    axes[0].set_xlim(ok.age_gyr.min() - 1, ok.age_gyr.max() + 1)
fig.suptitle(f"{LABEL}: the {len(ok)} fits inside 3σ ({ok.simulation_id.nunique()} runs); jittered ±0.03 Gyr",
             x=0.01, ha="left", fontsize=10.5, color=INK)
fig.savefig(f"{V}/analysis/passing_cells.png", dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"{V}/analysis/ks_vs_age.png\n{V}/analysis/passing_cells.png")
