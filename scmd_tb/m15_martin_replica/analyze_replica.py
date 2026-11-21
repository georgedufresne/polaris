#!/usr/bin/env python3
"""Collect the replica grid and summarise it: best KS by age, and, once null/ exists, the posterior age under the
weighting rules of Ying's papers (all on his KS statistic and his bootstrap null):
  - 2024/2025 2D KS, no cut, per isochrone ((run, age) unit): "probability" of the KS value in the null, read as KDE
    density, fitted-normal density or empirical survival fraction (the papers give no formula);
  - 2023 (M92) rule: one best fit per run, runs within 3 sigma of the null mean, fit probability 1 below the mean and
    S(x)/S(mean) above;
  - each times the calibration-star weight of his isochrones (../m15_ying_isochrones/ying_calib_chi2.csv, Eq. 5).
Writes <results>/../combined_<tag>.csv and prints to stdout (logs/analyze_replica*.out).
Options: --results DIR --null DIR --cal CSV --tag NAME (defaults: the first run: results/, null/, ying_calib_chi2.csv, "results")."""
import argparse
import glob

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, norm

E = "/dartfs-hpc/scratch/f007058/scmd_tb/m15_martin_replica"
ap = argparse.ArgumentParser()
ap.add_argument("--results", default=f"{E}/results")
ap.add_argument("--null", default=f"{E}/null")
ap.add_argument("--cal", default="/dartfs-hpc/scratch/f007058/scmd_tb/m15_ying_isochrones/ying_calib_chi2.csv")
ap.add_argument("--tag", default="results")
A = ap.parse_args()
cal = pd.read_csv(A.cal).set_index("simulation_id").p

g = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(f"{A.results}/task_*.csv"))], ignore_index=True)
g = g.drop_duplicates(["simulation_id", "age_myr"]).sort_values(["simulation_id", "age_myr"]).reset_index(drop=True)
g.to_csv(f"{E}/combined_{A.tag}.csv", index=False)
b = g.loc[g.ks_statistic.idxmin()]
print(f"{len(g)} units ({g.simulation_id.nunique()} runs); best KS {b.ks_statistic:.5f} at run {int(b.simulation_id)}, "
      f"{b.age_gyr:.1f} Gyr, DM {b.distance_modulus:.3f}, E {b.extinction_vi:.3f}")
print(f"fits at the box edges: E <= 0.081 {np.mean(g.extinction_vi <= 0.081):.0%}, DM <= 15.31 {np.mean(g.distance_modulus <= 15.31):.0%}, "
      f"DM >= 15.59 {np.mean(g.distance_modulus >= 15.59):.0%}")

try:
    null = np.sort(pd.concat([pd.read_csv(f) for f in glob.glob(f"{A.null}/null_task_*.csv")]).ks_statistic.values)
except ValueError:
    null = None
m, s = (null.mean(), null.std()) if null is not None else (np.nan, np.nan)
if null is not None:
    print(f"null: {len(null)} sets, mean {m:.5f}, sd {s:.5f}, max {null.max():.5f}, 3 sigma line {m + 3 * s:.5f}; best fit "
          f"{(b.ks_statistic - m) / s:.1f} sigma above the mean")

t = g.loc[g.groupby("age_gyr").ks_statistic.idxmin()].set_index("age_gyr")
print("\nbest fit at each age (his KS):  age  KS  [sigma above null]  run  DM  E   | median KS over runs")
med = g.groupby("age_gyr").ks_statistic.median()
for a, r in t.iterrows():
    sig = f"[{(r.ks_statistic - m) / s:5.1f}]" if null is not None else ""
    print(f"  {a:4.1f}  {r.ks_statistic:.5f} {sig}  {int(r.simulation_id)}  {r.distance_modulus:.3f}  {r.extinction_vi:.3f}   | {med[a]:.5f}")
if null is None:
    raise SystemExit("\nno null yet: posteriors skipped")


def wq(x, w, q):
    o = np.argsort(x)
    c = np.cumsum(w[o]) / w.sum()
    return np.interp(q, c, x[o])


def summary(name, d, logw):
    logw = np.asarray(logw, float)
    if not np.isfinite(logw).any():
        print(f"    {name:36s} every weight is zero")
        return
    w = np.exp(logw - logw[np.isfinite(logw)].max())
    a = d.age_gyr.values
    r = pd.Series(w).groupby(d.simulation_id.values).sum()
    print(f"    {name:36s} age {np.average(a, weights=w):.2f} +/- {np.sqrt(np.average((a - np.average(a, weights=w)) ** 2, weights=w)):.2f} "
          f"(median {wq(a, w, .5):.2f}, 16-84% {wq(a, w, .16):.2f}-{wq(a, w, .84):.2f}); DM {np.average(d.distance_modulus, weights=w):.3f}, "
          f"E {np.average(d.extinction_vi, weights=w):.3f}; effective runs {r.sum() ** 2 / (r ** 2).sum():.1f}")


kde = gaussian_kde(null)
S = lambda x: 1.0 - np.searchsorted(null, x, side="left") / len(null)     # noqa: E731
print("\nposterior age (Ying 2025, 2D KS: 13.39 +/- 0.39 Gyr):")
for label, d in (("per isochrone, (run, age) units, no cut", g),
                 ("per run, one best fit, no cut", g.loc[g.groupby("simulation_id").ks_statistic.idxmin()].reset_index(drop=True))):
    x = d.ks_statistic.values
    with np.errstate(divide="ignore"):
        lc = np.log(d.simulation_id.map(cal).values)
        rules = {"kde": kde.logpdf(x), "normal": norm.logpdf(x, m, s), "survival": np.log(S(x))}
    print(f"  {label}:")
    for k, lw in rules.items():
        summary(k, d, lw)
        summary(k + " x calibration", d, lw + lc)
d = g.loc[g.groupby("simulation_id").ks_statistic.idxmin()].reset_index(drop=True)
x = d.ks_statistic.values
inside = x <= m + 3 * s
with np.errstate(divide="ignore"):
    lw = np.log(np.where(x <= m, 1.0, S(x) / S(m)) * inside)
    lc = np.log(d.simulation_id.map(cal).values)
print(f"  M92 rule (per run, 3 sigma cut; {inside.sum()} of {len(d)} runs inside):")
summary("fit probability", d, lw)
summary("fit probability x calibration", d, lw + lc)
