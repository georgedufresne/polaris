#!/usr/bin/env python3
"""Are the DSEP4 47 Tuc isochrones usable for the sCMD fit?  One row per (run, isochrone file, age).

For each age 8.0-14.0 Gyr in isochrones_pyeep.txt (raw) and isochrones_example.txt (after
patch_lower_ms.py): number of points, lowest and highest initial mass, highest log L
(reaches the RGB if > 2), whether log L rises monotonically with mass up to the turnoff
region, and the largest gap in initial mass below 0.85 Msun (holes in the lower MS).
The sCMD needs F606W ~ 15-20 at mu ~ 13.4, i.e. masses ~0.55 Msun up the RGB.

    health.py FIRST LAST     -> out/health_FIRST_LAST.csv
"""
import os
import sys

import numpy as np
import pandas as pd

O = "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso_d4"
F = os.path.dirname(os.path.abspath(__file__))
AGES = np.round(np.arange(8.0, 14.01, 0.1), 1)


def read(path):
    rows = [l.split() for l in open(path).read().split("\n")[11:] if l.strip() and not l.lstrip().startswith("#")]
    a = np.array(rows, dtype=float)
    return pd.DataFrame({"age": a[:, 1] / 1e9, "m": a[:, 2], "logL": a[:, 7], "logT": a[:, 10]})


first, last = int(sys.argv[1]), int(sys.argv[2])
out = []
for r in range(first, last + 1):
    for kind in ("pyeep", "example"):
        p = f"{O}/{r}/isos/isochrones_{kind}.txt"
        if not os.path.isfile(p):
            out.append(dict(run=r, kind=kind, age=np.nan, status="missing"))
            continue
        try:
            d = read(p)
        except Exception as e:                       # unreadable file
            out.append(dict(run=r, kind=kind, age=np.nan, status=f"unreadable: {e}"[:80]))
            continue
        for a in AGES:
            b = d[np.isclose(d.age, a, atol=1e-3)].sort_values("m")
            if len(b) < 3:
                out.append(dict(run=r, kind=kind, age=a, status="no block", n=len(b)))
                continue
            low = b[b.m < 0.85]
            gap = np.diff(low.m.values).max() if len(low) > 1 else np.nan
            ms = b[b.logL < 0.2]                       # unevolved part: L should rise with mass
            mono = bool(np.all(np.diff(ms.logL.values) > -0.02)) if len(ms) > 2 else False
            out.append(dict(run=r, kind=kind, age=a, status="ok", n=len(b), m_min=b.m.min(), m_max=b.m.max(),
                            logL_max=b.logL.max(), logL_min=b.logL.min(), gap_low=gap, ms_monotonic=mono))
pd.DataFrame(out).to_csv(f"{F}/out/health_{first}_{last}.csv", index=False)
print(f"done {first}-{last}")
