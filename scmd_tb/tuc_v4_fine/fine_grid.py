#!/usr/bin/env python3
"""Does the v4 grid miss lower 2D KS values?  Two ways it could:

  1. stage 2 keeps only the 8 (DM, E) cells per (run, age) with the lowest KDE JS
     distance, and only those get a KS; the KS-best cell is often JS rank 7-8;
  2. the grid steps (DM 0.01, E(V-I) 0.005) may be coarse next to how fast the KS
     changes (one E step moves it 0.016 -> 0.05 at run 10091, 12.0 Gyr).

M. Ying minimised the KS itself over continuous (DM, reddening) (gp_minimize), so
neither applies to his fits.  For one (run, age) this scores the KS on

  coarse  every cell of the v4 grid (31 x 11), no JS pre-selection
  fine    DM +-0.02 in 0.0025, E +-0.005 in 0.0005 around the coarse minimum

with the v4 setup and the same common random numbers (RandomState(crn_seed) per cell).

    fine_grid.py RUN AGE [--quick]     (--quick: two cells only, to time and check replay)
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

F = os.path.dirname(os.path.abspath(__file__))
G = "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_v4"
sys.path.insert(0, os.path.join(G, "code"))
from isochronetoolbox.scmd.pipeline import Setup, context  # noqa: E402

run, age = int(sys.argv[1]), float(sys.argv[2])
s = Setup(**json.load(open(os.path.join(G, "setup.json"))))
ctx = context(s)
age_yr, make = ctx.model(run, age)


def ks(dm, e):
    cat = make(np.random.RandomState(s.crn_seed), dm, e)
    if cat is None:
        return np.nan
    return ctx.stage3(cat)["ks_statistic"]


if "--quick" in sys.argv:
    for dm, e in ((13.39, 0.03), (13.40, 0.03)):
        t = time.time()
        print(f"run {run} {age} Gyr  DM {dm}  E {e}:  KS {ks(dm, e):.5f}  ({time.time() - t:.1f} s)", flush=True)
    sys.exit()

rows = []
t0 = time.time()
for dm in s.dm:
    for e in s.ebv:
        rows.append(("coarse", round(dm, 4), round(e, 4), ks(dm, e)))
c = pd.DataFrame(rows, columns=["stage", "dm", "ebv", "ks"])
b = c.loc[c.ks.idxmin()]
print(f"coarse done in {time.time() - t0:.0f} s: min KS {b.ks:.5f} at DM {b.dm} E {b.ebv}", flush=True)
for dm in np.arange(b.dm - 0.02, b.dm + 0.02001, 0.0025):
    for e in np.arange(max(b.ebv - 0.005, 0.0), b.ebv + 0.00501, 0.0005):
        rows.append(("fine", round(dm, 4), round(e, 4), ks(dm, e)))
out = pd.DataFrame(rows, columns=["stage", "dm", "ebv", "ks"])
out.insert(0, "age_gyr", age)
out.insert(0, "simulation_id", run)
out.to_csv(os.path.join(F, "out", f"fine_{run}_{age:.1f}.csv"), index=False)
f = out[out.stage == "fine"].ks.idxmin()
print(f"fine done in {time.time() - t0:.0f} s: min KS {out.ks[f]:.5f} at DM {out.dm[f]} E {out.ebv[f]}")
