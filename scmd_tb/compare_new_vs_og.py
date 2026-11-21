#!/usr/bin/env python3
"""New toolbox sCMD pipeline vs the original (OG) results, same 500 MC runs.

Same estimator as analyze_scmd_only.py: rows with KS <= 0.0468, weight = KDE of
the bootstrap KS null at the row's KS, best (DM, E(B-V)) per (run, age),
weights normalised.  Reported for all 500 runs and for the 446 runs whose OG
stage 2 kept exactly four sCMDs per age (in the other 54 OG runs the stage-2
deletion did not happen for some ages, so stage 3 scored up to 308 per age).

Uncertainty on the new-minus-OG difference: cluster bootstrap over MC runs
(resample runs with replacement, recompute both weighted means).
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, ks_2samp

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
OUT = f"{B}/scmd_toolbox_vs_og_500runs.txt"
KS3 = 0.0468
boot_null = pd.read_csv(f"{B}/analysis_og500/Bootstrap/bootstrap_2d_ks_results.csv")["ks_statistic"].values
kde = gaussian_kde(boot_null, bw_method="scott")
ivars = pd.read_csv(f"{B}/analysis_og500/initial_vars.csv")
ivars["file_id"] = ivars["file_id"].astype(int)

og = pd.read_csv(f"{B}/analysis_og500/combined_2d_ks_results.csv")
new = pd.read_csv(f"{B}/og500/combined_2d_ks_results.csv")
runs = sorted(set(new.simulation_id) & set(og.simulation_id))
clean = [r for r, n in og.groupby("simulation_id").size().items() if n == 324 and r in runs]


def good(df):
    """The legacy prefilter's cut: KS read as float32, compared in float32."""
    return df.ks_statistic.astype(np.float32) <= KS3


def posterior_rows(df):
    f = df[good(df)].copy()
    f = f.loc[f.groupby(["simulation_id", "age_gyr"])["ks_statistic"].idxmin()]
    f["w"] = kde(f.ks_statistic.values)
    return f.merge(ivars, left_on="simulation_id", right_on="file_id", how="left")


def wstats(x, w):
    w = w / w.sum()
    m = np.sum(w * x)
    o = np.argsort(x)
    c = np.cumsum(w[o])
    q = lambda p: x[o][min(np.searchsorted(c, p), len(x) - 1)]
    return m, np.sqrt(np.sum(w * (x - m) ** 2)), q(0.16), q(0.5), q(0.84)


PARAMS = [("age_gyr", "Age (Gyr)"), ("FeH", "[Fe/H]"), ("distance_modulus", "DM"),
          ("extinction_ebv", "E(B-V)"), ("He Abundance", "Y"), ("CMIXLA", "alpha_MLT")]

lines = []
p = lines.append
p("sCMD pipeline: isochrone-toolbox (isochronetoolbox.scmd) vs the original three scripts")
p("=" * 92)
p("Same 500 MC runs (20000-20566), same OG isochrones (Isochrones_Only_July_10), same inputs;")
p("the toolbox reproduces the old scripts bit-for-bit for a given seed (tests/scmd/compare_with_legacy.py),")
p("so any difference here is Monte Carlo noise: the OG runs were unseeded.")
p("")
for label, sel in (("all 500 runs", runs), (f"{len(clean)} runs with a clean OG stage 2 (324 rows)", clean)):
    o = og[og.simulation_id.isin(sel)]
    n = new[new.simulation_id.isin(sel)]
    po, pn = posterior_rows(o), posterior_rows(n)
    p(f"--- {label} ---")
    p(f"  rows:              OG {len(o):>7,}   new {len(n):>7,}")
    p(f"  rows KS<=0.0468:   OG {good(o).sum():>7,}   new {good(n).sum():>7,}")
    p(f"  runs with a fit:   OG {po.simulation_id.nunique():>7}   new {pn.simulation_id.nunique():>7}   (of {len(sel)})")
    bo, bn = o.groupby("simulation_id").ks_statistic.min(), n.groupby("simulation_id").ks_statistic.min()
    both = bo.index.intersection(bn.index)
    d = (bn[both] - bo[both]).values
    p(f"  best KS per run:   OG median {bo.median():.4f}   new median {bn.median():.4f}   "
      f"paired diff {d.mean():+.4f} +/- {d.std(ddof=1) / np.sqrt(len(d)):.4f}   "
      f"2-sample KS p = {ks_2samp(bo.values, bn.values).pvalue:.2f}")
    p("")
    p(f"  {'parameter':<12} {'OG mean +/- std':>22} {'new mean +/- std':>22} {'new - OG (cluster bootstrap)':>34}")
    rs = np.random.RandomState(1)
    ids = np.array(sel)
    go, gn = dict(tuple(po.groupby("simulation_id"))), dict(tuple(pn.groupby("simulation_id")))
    for col, name in PARAMS:
        mo, so, *_ = wstats(po[col].values.astype(float), po.w.values)
        mn, sn, *_ = wstats(pn[col].values.astype(float), pn.w.values)
        diffs = []
        for _ in range(500):
            pick = ids[rs.randint(0, len(ids), len(ids))]
            a = pd.concat([go[r] for r in pick if r in go])
            b = pd.concat([gn[r] for r in pick if r in gn])
            diffs.append(np.average(b[col], weights=b.w) - np.average(a[col], weights=a.w))
        sd = np.std(diffs, ddof=1)
        p(f"  {name:<12} {mo:>12.3f} +/- {so:<7.3f} {mn:>12.3f} +/- {sn:<7.3f} "
          f"{mn - mo:>+12.3f} +/- {sd:.3f}  (z = {(mn - mo) / sd:+.2f})")
    p("")
open(OUT, "w").write("\n".join(lines) + "\n")
print("\n".join(lines))
print("saved", OUT)
