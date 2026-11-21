#!/usr/bin/env python3
"""47 Tuc: every statistic of the NGC 6791 method on the finished grid, plus two labelled variants.

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

sys.path.insert(0, "/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox")
from isochronetoolbox.scmd import posterior as P          # noqa: E402

G = sys.argv[1]
rows = pd.read_csv(f"{G}/combined_2d_ks_results.csv")
boot = pd.read_csv(f"{G}/bootstrap_2d_ks_results.csv")["ks_statistic"].values
mc = P.read_mc_inputs("/dartfs-hpc/scratch/f007058/47Tuc_var/varfeh-077.{run}", sorted(rows.simulation_id.unique()))
bins = P.BINS["47tuc"]()
mu, sd = boot.mean(), boot.std()
best = rows.ks_statistic.min()

import json
cell = json.load(open(f"{G}/bootstrap/cell.json"))
print(f"bootstrap null ({len(boot):,} realisations of run {cell['simulation_id']}, {cell['age_gyr']} Gyr, "
      f"mu {cell['distance_modulus']}, E(V-I) {cell['extinction_ebv']}):")
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
    d = f"{G}/posterior_{name}"
    os.makedirs(d, exist_ok=True)
    tab.to_csv(f"{d}/posterior_summary.csv", index=False)
    w.merge(mc, on="simulation_id", how="left").to_csv(f"{d}/posterior_weights.csv", index=False)
    for k, (c, h, _) in post.items():
        pd.DataFrame({"center": c, "probability": h}).to_csv(f"{d}/hist_{k}.csv", index=False)
    P.plot(post, f"{d}/posteriors.png", f"47 Tuc sCMD posteriors, {name} null: KS <= {lim:.4f}, "
           f"{len(w)} cells, {w.simulation_id.nunique()} runs", "E(V-I)")
    print(tab.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    top = w.sort_values("P_KS_normalized", ascending=False)
    print("  weight held by the top 1 / 5 / 10 runs: " + " / ".join(
        f"{top.groupby('simulation_id').P_KS_normalized.sum().sort_values(ascending=False).head(n).sum():.2f}" for n in (1, 5, 10)))
