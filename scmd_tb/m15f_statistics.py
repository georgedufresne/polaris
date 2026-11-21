#!/usr/bin/env python3
"""M15 (copy of tuc_statistics.py): every statistic of the NGC 6791 method on the finished grid, plus two labelled variants.

  strict   the NGC 6791 rule: KS <= bootstrap mean + 3 sigma, weights = null KDE at KS
  shifted  null moved so its mean is the best fit's KS: fits within 3 sampling-sigma of the best fit,
           weighted by how far above it they are (model-systematic floor, sampling width)
  scaled   null multiplied by best KS / null mean (width grows with the floor)

Writes posterior_<variant>/ (summary, weights, histograms, figure) in the grid directory and prints
the tables.
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/m15_grid_d4/code")   # the grid's own copy (identical to the home one): no Kerberos ticket needed
from isochronetoolbox.scmd import posterior as P          # noqa: E402

G = sys.argv[1]
rows = pd.read_csv(f"{G}/combined_2d_ks_results.csv")
EVI_MIN = float(os.environ.get("EVI_MIN", "-1"))          # e.g. 0.08: only cells at or above this E
SUFFIX = os.environ.get("POST_SUFFIX", "")
rows = rows[rows.extinction_ebv >= EVI_MIN - 1e-9]
if EVI_MIN > 0:
    print(f"cells with E >= {EVI_MIN} only: {len(rows)} rows")
# optional upper E and mu range (prior box), applied before the per-(run, age) minimum is taken
BOX = {k: float(os.environ[k]) for k in ("EVI_MAX", "MU_MIN", "MU_MAX") if k in os.environ}
if BOX:
    rows = rows[(rows.extinction_ebv <= BOX.get("EVI_MAX", 99) + 1e-9)
                & (rows.distance_modulus >= BOX.get("MU_MIN", -99) - 1e-9) & (rows.distance_modulus <= BOX.get("MU_MAX", 99) + 1e-9)]
    print(f"prior box {BOX} (and E >= {EVI_MIN}): {len(rows)} rows")
boot = pd.read_csv(f"{G}/bootstrap_2d_ks_results.csv")["ks_statistic"].values
mc = P.read_mc_inputs("/dartfs-hpc/scratch/f007058/M15/dsep4/var/varfeh-227.{run}", sorted(rows.simulation_id.unique()))
bins = P.BINS["47tuc"]()
# M15 final grid steps: 0.1 Gyr to 15 Gyr, mu 15.30-15.60 in 0.01, E 0.050-0.120 in 0.005
# (refined: E 0.000-0.120 in 0.001; bins stay 0.005 wide, from 0)
for k, (a, b, st) in {"age": (7.95, 15.15, 0.1), "feh": (-2.605, -1.895, 0.02), "dm": (15.295, 15.615, 0.01),
                      "ebv": (-0.0025, 0.1251, 0.005)}.items():
    e = np.arange(a, b, st)
    bins[k] = {"edges": e, "centers": (e[:-1] + e[1:]) / 2}
mu, sd = boot.mean(), boot.std()
best = rows.ks_statistic.min()

import json
cell = json.load(open(f"{G}/bootstrap/cell.json"))
print(f"bootstrap null ({len(boot):,} realisations of run {cell['simulation_id']}, {cell['age_gyr']} Gyr, "
      f"mu {cell['distance_modulus']}, E(606-814) {cell['extinction_ebv']}):")
print(f"  mean {mu:.5f}  std {sd:.5f}  median {np.median(boot):.5f}  99.9% {np.percentile(boot, 99.9):.5f}  max {boot.max():.5f}  3 sigma {mu + 3 * sd:.5f}")
best_run = rows.loc[rows.groupby("simulation_id").ks_statistic.idxmin()]
p_emp = [(boot >= k).mean() for k in best_run.ks_statistic]
print(f"  best fit KS {best:.5f} = {(best - mu) / sd:.1f} sigma above the null mean; empirical p of every run's best fit: "
      f"{'all < 1e-4' if max(p_emp) == 0 else max(p_emp)}")
print(f"  rows with KS <= 3 sigma: {(rows.ks_statistic <= mu + 3 * sd).sum()} of {len(rows)}")

variants = {"strict": boot, "shifted": boot + (best - mu), "scaled": boot * (best / mu)}
out = []
for name, null in variants.items():
    lim = null.mean() + 3 * null.std()
    w = P.weigh(rows, null, lim)
    print(f"\n[{name}] cut KS <= {lim:.5f}: {(rows.ks_statistic <= lim).sum()} rows, "
          f"{len(w)} (run, age) cells from {w.simulation_id.nunique() if len(w) else 0} runs carry weight")
    if w.empty or not np.isfinite(w.P_KS_normalized).all():
        print("  no posterior")
        continue
    post = P.posteriors(w, mc, bins)
    tab = P.summary_table(post, w)
    d = f"{G}/posterior_{name}{SUFFIX}"
    os.makedirs(d, exist_ok=True)
    tab.to_csv(f"{d}/posterior_summary.csv", index=False)
    w.merge(mc, on="simulation_id", how="left").to_csv(f"{d}/posterior_weights.csv", index=False)
    for k, (c, h, _) in post.items():
        pd.DataFrame({"center": c, "probability": h}).to_csv(f"{d}/hist_{k}.csv", index=False)
    P.plot(post, f"{d}/posteriors.png", f"M15 sCMD posteriors, {name} null: KS <= {lim:.4f}, "
           f"{len(w)} cells, {w.simulation_id.nunique()} runs", "E(F606W-F814W) = 1.02 E(B-V)")
    print(tab.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    top = w.sort_values("P_KS_normalized", ascending=False)
    print("  weight held by the top 1 / 5 / 10 runs: " + " / ".join(
        f"{top.groupby('simulation_id').P_KS_normalized.sum().sort_values(ascending=False).head(n).sum():.2f}" for n in (1, 5, 10)))
