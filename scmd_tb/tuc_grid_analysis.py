#!/usr/bin/env python3
"""47 Tuc test grid: what fits, and with which physics.

Reads the collected grid (combined_2d_ks_results.csv, kde_js_scores.csv.gz),
drops ages a run's tracks cannot support (tuc_age_limits.json), and writes to
tuc_grid/analysis/:

  best_per_run.csv          best (age, mu, E(V-I)) of every run by 2D KS, with its MC inputs
  top_runs/rankNN_*.png     observed CMD next to the best sCMD of the 10 best runs
  top10_sheet.png           the ten best sCMDs on one page
  ks_vs_age.png             min KS over (mu, E(V-I)) at each age, every run
  best_fit_distributions.png  best-fit age, mu, E(V-I) over the runs
  ks_vs_physics.png         best KS of each run against its Monte Carlo inputs
  js_map_best_run.png       stage-2 JS over (age, mu) for the best run, minimised over E(V-I)
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
G = sys.argv[1] if len(sys.argv) > 1 else f"{B}/tuc_grid"
# replay with the grid's own frozen code, so regenerated sCMDs are the ones scored
sys.path.insert(0, f"{G}/code" if os.path.isdir(f"{G}/code/isochronetoolbox") else "/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox")
from isochronetoolbox.scmd import context, preset, regenerate          # noqa: E402
from isochronetoolbox.scmd.plot import cmd_pair                          # noqa: E402
A = f"{G}/analysis"
os.makedirs(f"{A}/top_runs", exist_ok=True)

lim = {int(k): v for k, v in json.load(open(f"{B}/tuc_age_limits.json")).items()}
ks = pd.read_csv(f"{G}/combined_2d_ks_results.csv")
js = pd.read_csv(f"{G}/kde_js_scores.csv.gz")
ks["valid"] = ks.age_gyr < ks.simulation_id.map(lim) - 1e-6
js["valid"] = js.age_gyr < js.simulation_id.map(lim) - 1e-6
print(f"{ks.simulation_id.nunique()} runs, {len(ks):,} stage-3 rows ({(~ks.valid).sum():,} at ages past a run's track coverage, dropped)")
ks, js = ks[ks.valid], js[js.valid]

# Monte Carlo inputs of each run
mc = []
for r in sorted(ks.simulation_id.unique()):
    v = dict((l.split()[1], float(l.split()[0])) for l in open(f"/dartfs-hpc/scratch/f007058/47Tuc_var/varfeh-077.{r}") if l.strip())
    v["simulation_id"] = r
    mc.append(v)
mc = pd.DataFrame(mc)

best = ks.loc[ks.groupby("simulation_id").ks_statistic.idxmin()].sort_values("ks_statistic")
best = best.merge(mc, on="simulation_id")
best.to_csv(f"{A}/best_per_run.csv", index=False)
cols = ["simulation_id", "age_gyr", "distance_modulus", "extinction_ebv", "ks_statistic", "js_divergence",
        "FeH", "Yprim", "CMIXLA", "FGRY", "ALPHAE", "ALPHAC"]
print(best[cols].head(15).rename(columns={"distance_modulus": "mu", "extinction_ebv": "E(V-I)"}).to_string(index=False))

# ---- observed vs sCMD for the ten best runs --------------------------------------------
import json as _json
from isochronetoolbox.scmd import Setup
_saved = _json.load(open(f"{G}/setup.json"))
s = Setup(**{k: v for k, v in _saved.items() if k in Setup.__dataclass_fields__})
ctx = context(s)
(c0, c1), (m0, m1) = s.ks_window
top = best.head(10)
cats = []
for k, r in enumerate(top.itertuples(), 1):
    cat = regenerate(ctx, int(r.simulation_id), r.age_gyr, r.distance_modulus, r.extinction_ebv)
    chk = ctx.stage3(cat)["ks_statistic"]
    note = "" if abs(chk - r.ks_statistic) < 1e-9 else f"  (replay KS {chk:.4f})"
    title = (f"#{k}  47 Tuc  run {int(r.simulation_id)}   {r.age_gyr:.1f} Gyr   mu {r.distance_modulus:.2f}   "
             f"E(V-I) {r.extinction_ebv:.3f}\nKS {r.ks_statistic:.4f}   JS {r.js_divergence:.4f}   "
             f"[Fe/H] {r.FeH:.3f}  Y {r.Yprim:.3f}  alpha_MLT {r.CMIXLA:.3f}{note}")
    cmd_pair(ctx.obs, ctx.points(cat), title, f"{A}/top_runs/rank{k:02d}_run{int(r.simulation_id)}.png",
             window=s.ks_window, labels=("F606W - F814W", "F606W"), xlim=(0.38, 0.86), ylim=(m1 + 0.15, m0 - 0.15))
    cats.append((r, ctx.points(cat)))

fig, axes = plt.subplots(2, 6, figsize=(22, 8.5), sharex=True, sharey=True)
axes = axes.ravel()
axes[0].scatter(ctx.obs[:, 0], ctx.obs[:, 1], s=0.15, c="k", alpha=0.3, linewidths=0, rasterized=True)
axes[0].set_title(f"Observed (N = {len(ctx.obs)})", fontsize=10)
for ax, (r, p) in zip(axes[1:], cats):
    ax.scatter(p[:, 0], p[:, 1], s=0.15, c="#c0392b", alpha=0.3, linewidths=0, rasterized=True)
    ax.set_title(f"run {int(r.simulation_id)}  {r.age_gyr:.1f} Gyr  KS {r.ks_statistic:.3f}", fontsize=10)
axes[-1].axis("off")
for ax in axes:
    ax.set_xlim(0.38, 0.86)
    ax.set_ylim(m1 + 0.15, m0 - 0.15)
fig.supxlabel("F606W - F814W")
fig.supylabel("F606W")
fig.suptitle("47 Tuc: observed CMD and the best sCMD of the ten best-fitting Monte Carlo runs", fontsize=13)
fig.tight_layout()
fig.savefig(f"{A}/top10_sheet.png", dpi=130)
plt.close(fig)

# ---- KS vs age ----------------------------------------------------------------------------
per_age = ks.groupby(["simulation_id", "age_gyr"]).ks_statistic.min().reset_index()
fig, ax = plt.subplots(figsize=(9, 5.5))
for r, g in per_age.groupby("simulation_id"):
    hl = r in set(top.simulation_id)
    ax.plot(g.age_gyr, g.ks_statistic, color="#c0392b" if hl else "0.75", lw=1.4 if hl else 0.6, zorder=3 if hl else 1)
ax.set_xlabel("isochrone age (Gyr)")
ax.set_ylabel("min 2D KS over (mu, E(V-I))")
ax.set_title("47 Tuc: KS against age for every run (red: the ten best runs)")
ax.grid(alpha=0.3)
fig.tight_layout()
fig.savefig(f"{A}/ks_vs_age.png", dpi=140)
plt.close(fig)

# ---- best-fit distributions ----------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
for ax, (c, lab, bins) in zip(axes, [("age_gyr", "best-fit age (Gyr)", np.arange(7.95, 14.1, 0.2)),
                                      ("distance_modulus", "best-fit mu", np.arange(13.145, 13.46, 0.01)),
                                      ("extinction_ebv", "best-fit E(V-I)", np.arange(-0.0025, 0.055, 0.005))]):
    ax.hist(best[c], bins=bins, color="0.6", label="all runs")
    ax.hist(top[c], bins=bins, color="#c0392b", label="ten best")
    ax.set_xlabel(lab)
axes[0].legend()
fig.suptitle("47 Tuc: best fit of each run")
fig.tight_layout()
fig.savefig(f"{A}/best_fit_distributions.png", dpi=140)
plt.close(fig)

# ---- KS vs physics ------------------------------------------------------------------------
phys = [("FeH", "[Fe/H]"), ("Yprim", "Y"), ("CMIXLA", "alpha_MLT"), ("FGRY", "diffusion (He)"),
        ("ALPHAE", "envelope overshoot"), ("ALPHAC", "core overshoot"), ("alphafe", "[alpha/Fe]"), ("KTTAU", "T-tau")]
fig, axes = plt.subplots(2, 4, figsize=(16, 7.5), sharey=True)
for ax, (c, lab) in zip(axes.ravel(), phys):
    ax.scatter(best[c], best.ks_statistic, s=14, c="0.4")
    ax.scatter(top[c], top.ks_statistic, s=22, c="#c0392b")
    rho = pd.Series(best[c]).corr(best.ks_statistic, method="spearman")
    ax.set_title(f"{lab}   (Spearman {rho:+.2f})", fontsize=10)
axes[0, 0].set_ylabel("best KS of the run")
axes[1, 0].set_ylabel("best KS of the run")
fig.suptitle("47 Tuc: how well each run fits against its Monte Carlo inputs (red: ten best)")
fig.tight_layout()
fig.savefig(f"{A}/ks_vs_physics.png", dpi=140)
plt.close(fig)

# ---- JS map for the best run --------------------------------------------------------------
r0 = int(top.simulation_id.iloc[0])
j = js[js.simulation_id == r0].groupby(["age_gyr", "distance_modulus"]).js_distance.min().unstack()
fig, ax = plt.subplots(figsize=(8, 5.5))
im = ax.pcolormesh(j.columns, j.index, j.values, shading="nearest", cmap="viridis_r")
fig.colorbar(im, label="stage-2 KDE JS distance (min over E(V-I))")
ax.set_xlabel("mu")
ax.set_ylabel("age (Gyr)")
ax.set_title(f"47 Tuc run {r0}: the age-distance valley")
fig.tight_layout()
fig.savefig(f"{A}/js_map_best_run.png", dpi=140)
plt.close(fig)
print("figures in", A)
