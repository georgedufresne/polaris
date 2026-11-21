#!/usr/bin/env python3
"""Ying et al. (2024, 2025) 2D-KS weighting with no cut, on the finished KS tables of this run and of
../m15_yingiso_nozp.

The papers say only that each isochrone gets "the probability" of its KS value in the bootstrap reference
distribution (Ying et al. 2024 Sec. 5.4, "p(L=L*|theta-hat)"; 2025 Sec. IV.2 and Fig. 8, null shown with a KDE), with no
threshold and no formula. Three readings, none cut:
  kde      density of a Gaussian KDE of the null at the KS value (scipy, Scott bandwidth; log density, so far tails
           do not underflow)
  normal   density of a normal fitted to the null (mean, std)
  survival empirical fraction of null realisations >= the KS value (zero beyond the null's maximum)
Each is applied per isochrone, i.e. per (run, age) cell (best over mu, E), as "each isochrone tested" in the 2025
paper, and per run (one best fit per run), with and without calibration-star weights for his isochrones
(../m15_ying_isochrones/ying_calib_chi2.csv). Output: logs/nocut_posterior.out."""
import json

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, norm

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
LIM = {int(k): v for k, v in json.load(open("/dartfs-hpc/scratch/f007058/M15/iso_ying2025/age_limits_ying.json")).items()}
cal = pd.read_csv(f"{B}/m15_ying_isochrones/ying_calib_chi2.csv").set_index("simulation_id").p


def wq(x, w, q):
    o = np.argsort(x)
    c = np.cumsum(w[o]) / w.sum()
    return np.interp(q, c, x[o])


def summary(name, d, logw):
    logw = np.asarray(logw, float)
    if not np.isfinite(logw).any():
        print(f"    {name:34s} every weight is zero")
        return
    w = np.exp(logw - logw[np.isfinite(logw)].max())
    a = d.age_gyr.values
    r = pd.Series(w).groupby(d.simulation_id.values).sum()
    top = d.iloc[np.argsort(-w)[:3]]
    print(f"    {name:34s} age {np.average(a, weights=w):.2f} (median {wq(a, w, .5):.2f}, 16-84% {wq(a, w, .16):.2f}-"
          f"{wq(a, w, .84):.2f}); DM {np.average(d.distance_modulus, weights=w):.3f}, E {np.average(d.extinction_ebv, weights=w):.3f}; "
          f"effective runs {r.sum() ** 2 / (r ** 2).sum():.1f}; weight in top cell {w.max() / w.sum():.0%}; top ages "
          + ", ".join(f"{t.simulation_id}:{t.age_gyr:.1f}" for t in top.itertuples()))


for run in ("m15_yingiso_nozp", "m15_yingiso_err99"):
    d = f"{B}/{run}_refined"
    g = pd.read_csv(f"{d}/combined_2d_ks_results.csv",
                    usecols=["simulation_id", "age_gyr", "distance_modulus", "extinction_ebv", "ks_statistic"])
    g = g[g.age_gyr < g.simulation_id.map(LIM) - 1e-6].reset_index(drop=True)
    boot = np.sort(pd.read_csv(f"{d}/bootstrap_2d_ks_results.csv").ks_statistic.values)
    m, s = boot.mean(), boot.std()
    kde = gaussian_kde(boot)
    print(f"\n{run}: null {m:.5f} +/- {s:.5f} (max {boot.max():.5f}, KDE bandwidth {kde.factor * s:.5f}); "
          f"best KS {g.ks_statistic.min():.5f} ({(g.ks_statistic.min() - m) / s:.1f} sigma)")
    for label, t in (("per isochrone, (run, age) cells", g.loc[g.groupby(["simulation_id", "age_gyr"]).ks_statistic.idxmin()]),
                     ("per run, one best fit", g.loc[g.groupby("simulation_id").ks_statistic.idxmin()])):
        t = t.reset_index(drop=True)
        x = t.ks_statistic.values
        with np.errstate(divide="ignore"):
            logs = {"kde": kde.logpdf(x), "normal": norm.logpdf(x, m, s),
                    "survival": np.log(1.0 - np.searchsorted(boot, x, side="left") / len(boot))}
            lc = np.log(t.simulation_id.map(cal).values)     # some calibration weights underflow to 0 -> weight 0
        print(f"  {label} ({len(t)}):")
        for k, lw in logs.items():
            summary(k, t, lw)
            summary(k + " x calibration", t, lw + lc)
