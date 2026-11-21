#!/usr/bin/env python3
"""Would M. Ying's version of the 2D KS statistic pass where ours fails?

Same sCMDs (v4 model D), four statistics of |ECDF_obs - ECDF_sim| at the observed points:

  ours       ECDF quadrant (bluer, brighter): rank2d on (colour, F606W); 95th percentile
  flip       his quadrant (bluer, fainter): rank2d on (colour, -F606W); 95th percentile
  p99mask    our quadrant; 99th percentile of the points where neither ECDF is 0
  martin     his quadrant, 99th percentile, zero points dropped (vorbin_chi2.KS_2d)

Each statistic gets its own null, built exactly as the v4 bootstrap (fresh reference
RandomState(2,000,000 + i), artificial stars resampled RandomState(1,000,000 + i),
sCMD RandomState(i)) at the v4 best cell, so the pass/fail cut is like for like.

    martin_stat.py fit                    the fits (observed stars vs the grid sCMD)
    martin_stat.py null FIRST N           realisations FIRST .. FIRST+N-1 of the null
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay

F = os.path.dirname(os.path.abspath(__file__))
G = "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_v4"
sys.path.insert(0, os.path.join(G, "code"))
from isochronetoolbox.scmd.pipeline import Setup, _window, context  # noqa: E402
from isochronetoolbox.scmd.stats import rank2d  # noqa: E402

s = Setup(**json.load(open(os.path.join(G, "setup.json"))))
ctx = context(s)
CELL = (10091, 12.0, 13.39, 0.03)                  # v4 best cell: where the null is built


def stats(ref, sim):
    """The four statistics of reference points ``ref`` against sCMD points ``sim``."""
    tri = Delaunay(sim)
    out = {}
    for flip in (False, True):
        f = np.array([1.0, -1.0 if flip else 1.0])
        r_cdf = rank2d(ref * f) / len(ref)
        s_cdf = LinearNDInterpolator(tri, rank2d(sim * f) / len(sim), fill_value=0)(ref)
        d = np.abs(r_cdf - s_cdf)
        m = (r_cdf != 0) & (s_cdf != 0)
        dm = np.sort(d[m])
        p99 = dm[int(len(dm) * 0.99)]                # np.partition(..., int(n * 0.99)) in vorbin_chi2
        if flip:
            out["flip"], out["martin"] = np.percentile(d, 95), p99
        else:
            out["ours"], out["p99mask"] = np.percentile(d, 95), p99
    return out


def points(cat):
    p = ctx.points(cat)
    return p[_window(p, s.ks_window)]


if sys.argv[1] == "fit":
    cells = [CELL, (10091, 12.0, 13.40, 0.03), (10046, 12.5, 13.32, 0.0), (10099, 10.9, 13.43, 0.02)]
    rows = []
    for run, age, dm, e in cells:
        _, make = ctx.model(run, age)
        r = stats(ctx.obs3, points(make(np.random.RandomState(s.crn_seed), dm, e)))
        rows.append(dict(run=run, age=age, dm=dm, ebv=e, **r))
        print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(os.path.join(F, "out", "martin_stat_fit.csv"), index=False)
else:
    first, n = int(sys.argv[2]), int(sys.argv[3])
    run, age, dm, e = CELL
    _, make = ctx.model(run, age)
    full = ctx.ast
    rows = []
    for i in range(first, first + n):
        ctx.ast = full
        ref = points(make(np.random.RandomState(2_000_000 + i), dm, e))
        ctx.ast = full.resample(1_000_000 + i)
        sim = points(make(np.random.RandomState(i), dm, e))
        rows.append(dict(i=i, **stats(ref, sim)))
    ctx.ast = full
    pd.DataFrame(rows).to_csv(os.path.join(F, "out", f"martin_stat_null_{first:05d}.csv"), index=False)
    print(f"done {first}..{first + n - 1}")
