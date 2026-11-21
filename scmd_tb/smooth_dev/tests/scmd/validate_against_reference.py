"""Validate a finished sCMD grid against reference results on the same runs.

The third tier of the sCMD checks (see the README):

  test_scmd.py              seconds     building blocks, CI
  compare_with_legacy.py    minutes     bit-for-bit against the old scripts, one (run, age) at a time
  this script               a grid      the science output: posteriors from two result sets agree

It applies the estimator of the original analysis (analyze_scmd_only.py) to
both sets: rows with KS <= 0.0468 (cut in float32, as the old prefilter did),
the best (DM, E(B-V)) per (run, age), each weighted by the density of the
bootstrap KS null at its KS.  The difference of the weighted means gets an
uncertainty from a cluster bootstrap over runs, and the check fails if any
|z| exceeds --max-z.

With the same isochrones and settings the two sets differ only by Monte Carlo
noise, so this checks that a new or changed pipeline still gives the same
science.  With different isochrones (e.g. emulated vs DSEP) it measures the
shift instead -- pass --max-z 0 to report without failing.

  python tests/scmd/validate_against_reference.py --new OUT \\
      --ref /dartfs-hpc/rc/lab/C/ChaboyerB/2d_ks [--runs runs.txt] [--report FILE]

``--new``/``--ref`` take an output directory of this package (its
combined_2d_ks_results.csv), a directory of legacy 2d_ks_results_<run>.csv
files, or one combined CSV.  500 runs took 417 CPU-hours (1 h 44 min on 272
CPUs); 50 runs give the age difference to about +/- 0.04 Gyr.
"""

from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data", "ngc6791")
KS3 = 0.0468
PARAMS = [("age_gyr", "Age (Gyr)"), ("distance_modulus", "DM"), ("extinction_ebv", "E(B-V)"),
          ("FeH", "[Fe/H]"), ("He Abundance", "Y"), ("CMIXLA", "alpha_MLT")]
COLS = ["simulation_id", "age_gyr", "distance_modulus", "extinction_ebv", "ks_statistic"]


def load(path: str, runs=None) -> pd.DataFrame:
    """Stage-3 rows from an output dir, a legacy 2d_ks dir, or a combined CSV."""
    if os.path.isdir(path):
        comb = os.path.join(path, "combined_2d_ks_results.csv")
        if os.path.exists(comb):
            df = pd.read_csv(comb, usecols=COLS)
        else:
            files = glob.glob(os.path.join(path, "2d_ks_results_*.csv"))
            if runs is not None:
                keep = {f"2d_ks_results_{r}.csv" for r in runs}
                files = [f for f in files if os.path.basename(f) in keep]
            df = pd.concat([pd.read_csv(f, usecols=COLS) for f in files], ignore_index=True)
    else:
        df = pd.read_csv(path, usecols=COLS)
    return df if runs is None else df[df.simulation_id.isin(runs)]


def posterior_rows(df: pd.DataFrame, null_kde, ivars: pd.DataFrame) -> pd.DataFrame:
    f = df[df.ks_statistic.astype(np.float32) <= KS3]
    f = f.loc[f.groupby(["simulation_id", "age_gyr"])["ks_statistic"].idxmin()].copy()
    f["w"] = null_kde(f.ks_statistic.values)
    return f.merge(ivars, left_on="simulation_id", right_on="file_id", how="left")


def per_run_sums(p: pd.DataFrame, runs, cols):
    """Sum of weights and of weight x parameter per run (zero for runs without fits)."""
    g = p.groupby("simulation_id")
    W = g["w"].sum().reindex(runs, fill_value=0.0).values
    S = np.column_stack([(p["w"] * p[c]).groupby(p.simulation_id).sum().reindex(runs, fill_value=0.0).values
                         for c in cols])
    return W, S


def wstats(x, w):
    w = w / w.sum()
    m = np.sum(w * x)
    o = np.argsort(x)
    c = np.cumsum(w[o])
    q = lambda f: x[o][min(np.searchsorted(c, f), len(x) - 1)]
    return m, np.sqrt(np.sum(w * (x - m) ** 2)), q(0.16), q(0.84)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--new", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--runs", help="file of run numbers (default: runs present in both)")
    ap.add_argument("--max-z", type=float, default=3.0, help="fail above this |z|; 0 = report only")
    ap.add_argument("--n-boot", type=int, default=2000)
    ap.add_argument("--report", help="also write the table to this file")
    a = ap.parse_args()

    runs = [int(x) for x in open(a.runs).read().split()] if a.runs else None
    new, ref = load(a.new, runs), load(a.ref, runs)
    runs = sorted(set(new.simulation_id) & set(ref.simulation_id))
    new, ref = new[new.simulation_id.isin(runs)], ref[ref.simulation_id.isin(runs)]
    null = gaussian_kde(pd.read_csv(os.path.join(DATA, "ks_bootstrap_null.csv"))["ks_statistic"].values,
                        bw_method="scott")
    ivars = pd.read_csv(os.path.join(DATA, "mc_initial_vars.csv"))
    ivars["file_id"] = ivars["file_id"].astype(int)
    pn, pr = posterior_rows(new, null, ivars), posterior_rows(ref, null, ivars)
    cols = [c for c, _ in PARAMS if c in pn.columns and pn[c].notna().all() and pr[c].notna().all()]

    Wn, Sn = per_run_sums(pn, runs, cols)
    Wr, Sr = per_run_sums(pr, runs, cols)
    rs = np.random.RandomState(12345)
    counts = rs.multinomial(len(runs), np.full(len(runs), 1.0 / len(runs)), size=a.n_boot)
    boot = (counts @ Sn) / (counts @ Wn)[:, None] - (counts @ Sr) / (counts @ Wr)[:, None]

    out = []
    p = out.append
    p(f"{len(runs)} runs in both sets")
    p(f"  rows:            new {len(new):>8,}   reference {len(ref):>8,}")
    p(f"  rows KS<={KS3}:  new {(new.ks_statistic.astype(np.float32) <= KS3).sum():>8,}"
      f"   reference {(ref.ks_statistic.astype(np.float32) <= KS3).sum():>8,}")
    p(f"  runs with a fit: new {pn.simulation_id.nunique():>8}   reference {pr.simulation_id.nunique():>8}")
    p("")
    p(f"  {'parameter':<11}{'new: mean +/- std [p16, p84]':>36}{'reference':>36}{'new - ref':>22}{'z':>7}")
    worst = 0.0
    for k, c in enumerate(cols):
        name = dict(PARAMS)[c]
        mn, sn, ln, hn = wstats(pn[c].values.astype(float), pn.w.values)
        mr, sr, lr, hr = wstats(pr[c].values.astype(float), pr.w.values)
        d, e = mn - mr, np.std(boot[:, k], ddof=1)
        z = d / e if e > 0 else 0.0
        worst = max(worst, abs(z))
        p(f"  {name:<11}{mn:>12.3f} +/- {sn:<6.3f}[{ln:.3f}, {hn:.3f}]"
          f"{mr:>12.3f} +/- {sr:<6.3f}[{lr:.3f}, {hr:.3f}]{d:>+12.3f} +/- {e:<6.3f}{z:>+7.2f}")
    ok = a.max_z <= 0 or worst <= a.max_z
    p("")
    p(f"largest |z| = {worst:.2f}" + ("" if a.max_z <= 0 else f"  ->  {'PASS' if ok else 'FAIL'} (limit {a.max_z})"))
    text = "\n".join(out)
    print(text)
    if a.report:
        open(a.report, "w").write(text + "\n")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
