#!/usr/bin/env python3
"""Oldest age each DSEP4 run's isochrones support, for the bootstrap cell and the figures
(M15 copy of scmd_tb/tuc_d4v2_build/age_limits_d4.py).

An age block of isochrones_example.txt is usable if it has >= 30 points, reaches the RGB
(max log L > 1.5) and starts on the lower main sequence (lowest initial mass < 0.6 Msun).
The limit is the first unusable age from 8.0 Gyr up (16.1 if all 8.0-16.0 are usable), so
ages below it pass ``age < limit``.  Prints how many runs are limited, and where.

    age_limits_d4.py RUNS_FILE OUT_JSON
"""
import json
import sys

import numpy as np
import pandas as pd

O = "/dartfs-hpc/scratch/f007058/M15/iso_dsep4"


def blocks(path):
    rows = [l.split() for l in open(path).read().split("\n")[11:] if l.strip() and not l.lstrip().startswith("#")]
    a = np.array(rows, dtype=float)
    d = pd.DataFrame({"age": np.round(a[:, 1] / 1e9, 1), "m": a[:, 2], "logL": a[:, 7]})
    return d.groupby("age").agg(n=("m", "size"), mmin=("m", "min"), lmax=("logL", "max"))


runs = [int(l) for l in open(sys.argv[1]) if l.strip()]
lim = {}
for r in runs:
    g = blocks(f"{O}/{r}/isos/isochrones_example.txt").reindex(np.round(np.arange(8.0, 16.01, 0.1), 1))
    bad = ~((g.n >= 30) & (g.lmax > 1.5) & (g.mmin < 0.6))
    lim[r] = float(g.index[bad.values][0]) if bad.any() else 16.1
json.dump({str(k): v for k, v in sorted(lim.items())}, open(sys.argv[2], "w"), indent=0)
v = np.array(list(lim.values()))
print(f"age limits for {len(v)} runs: {np.sum(v > 16)} usable to 16 Gyr; {np.sum(v <= 14.0)} limited at or below 14 Gyr "
      f"(the grid's top age); lowest {v.min():.1f}; written to {sys.argv[2]}")
