#!/usr/bin/env python3
"""Ying et al. (2023, M92) weighting rule applied to the finished KS tables of this run and of ../m15_yingiso_err99.

The M92 rule (Sec. 4.1-4.3 of Ying, Chaboyer, Boudreaux et al. 2023, AJ 166, 18):
  - one best fit per MC set: the minimum statistic over (age, mu, E) gives that set's age;
  - keep sets within 3 sigma of the mean of the null (bootstrap) distribution;
  - "fit probability" 1 below the null mean; above it, from the null's distribution function. The paper gives no
    formula, so two readings are reported: S(x)/S(mean) (continuous at the mean) and S(x) (as literally stated),
    with S the empirical survival fraction of the bootstrap null;
  - final weight = fit probability x calibration-star weight (per set).
For comparison: per-set best with our weights (null KDE density), and our strict posterior (every (run, age) cell
within 3 sigma, weights = null KDE). Calibration weights: ../m15_ying_isochrones/ying_calib_chi2.csv (his isochrones).
Output: logs/m92_rule_posterior.out."""
import json

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
LIM = {int(k): v for k, v in json.load(open("/dartfs-hpc/scratch/f007058/M15/iso_ying2025/age_limits_ying.json")).items()}
cal = pd.read_csv(f"{B}/m15_ying_isochrones/ying_calib_chi2.csv").set_index("simulation_id").p


def wq(x, w, q):
    o = np.argsort(x)
    c = np.cumsum(w[o]) / w.sum()
    return np.interp(q, c, x[o])


def summary(name, d, w):
    w = np.asarray(w, float)
    if w.sum() <= 0:
        print(f"  {name:44s} no weight (nothing within 3 sigma)")
        return
    a = d.age_gyr.values
    r = pd.Series(w).groupby(d.simulation_id.values).sum()
    print(f"  {name:44s} age {np.average(a, weights=w):.2f} (median {wq(a, w, .5):.2f}, 16-84% {wq(a, w, .16):.2f}-"
          f"{wq(a, w, .84):.2f}); DM {np.average(d.distance_modulus, weights=w):.3f}, E {np.average(d.extinction_ebv, weights=w):.3f}; "
          f"{(r > 0).sum()} runs, effective {r.sum() ** 2 / (r ** 2).sum():.0f}")


for run in ("m15_yingiso_nozp", "m15_yingiso_err99"):
    d = f"{B}/{run}_refined"
    g = pd.read_csv(f"{d}/combined_2d_ks_results.csv",
                    usecols=["simulation_id", "age_gyr", "distance_modulus", "extinction_ebv", "ks_statistic"])
    g = g[g.age_gyr < g.simulation_id.map(LIM) - 1e-6].reset_index(drop=True)
    boot = np.sort(pd.read_csv(f"{d}/bootstrap_2d_ks_results.csv").ks_statistic.values)
    m, s = boot.mean(), boot.std()
    S = lambda x: 1.0 - np.searchsorted(boot, x, side="right") / len(boot)     # noqa: E731
    best = g.loc[g.groupby("simulation_id").ks_statistic.idxmin()].reset_index(drop=True)
    inside = best.ks_statistic <= m + 3 * s
    print(f"\n{run}: null {m:.5f} +/- {s:.5f}, 3 sigma line {m + 3 * s:.5f}; {len(best)} runs, best KS {best.ks_statistic.min():.5f} "
          f"({(best.ks_statistic.min() - m) / s:.1f} sigma); runs whose best fit is within 3 sigma: {inside.sum()}")
    x = best.ks_statistic.values
    p_cont = np.where(x <= m, 1.0, S(x) / S(m)) * inside
    p_lit = np.where(x <= m, 1.0, S(x)) * inside
    p_kde = gaussian_kde(boot)(x) * inside
    c = best.simulation_id.map(cal).values
    print(" one best fit per run (M92 rule):")
    summary("fit probability S(x)/S(mean)", best, p_cont)
    summary("fit probability S(x) (literal)", best, p_lit)
    summary("our weights (null KDE density)", best, p_kde)
    summary("S(x)/S(mean) x calibration weights", best, p_cont * c)
    summary("S(x) x calibration weights", best, p_lit * c)
    cells = g.loc[g.groupby(["simulation_id", "age_gyr"]).ks_statistic.idxmin()].reset_index(drop=True)
    xi = cells.ks_statistic.values
    pin = xi <= m + 3 * s
    print(" every (run, age) cell (ours):")
    summary("strict: our weights (null KDE density)", cells, gaussian_kde(boot)(xi) * pin)
    summary("strict: S(x)/S(mean) weights", cells, np.where(xi <= m, 1.0, S(xi) / S(m)) * pin)
