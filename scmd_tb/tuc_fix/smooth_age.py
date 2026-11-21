#!/usr/bin/env python3
"""Best age per run from a smooth fit to the KS-vs-age curve, not the lowest grid point.

For each run: the lowest KS over (DM, E) at each age; a quadratic in age fitted to the
ages within +-HALF Gyr of the grid minimum (re-centred once on the fitted vertex);
the vertex is the smooth best age. HALF = 1.2 Gyr spans about four sawtooth periods.
"""
import sys
import numpy as np
import pandas as pd

F = "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_fix"
HALF = 1.2


def smooth_age(c):
    a, y = c.index.values, c.values
    centre = a[np.argmin(y)]
    for _ in range(2):
        s = np.abs(a - centre) <= HALF
        p = np.polyfit(a[s], y[s], 2)
        if p[0] <= 0:
            return np.nan, np.nan
        centre = np.clip(-p[1] / (2 * p[0]), a.min(), a.max())
    return centre, np.polyval(p, centre)


def table(v):
    g = pd.read_csv(f"{F}/{v}/combined_2d_ks_results.csv")
    out = {}
    for r, gr in g.groupby("simulation_id"):
        c = gr.groupby("age_gyr").ks_statistic.min()
        out[r] = (c.idxmin(), *smooth_age(c))
    return pd.DataFrame(out, index=["grid_age", "smooth_age", "smooth_ks"]).T


V = sys.argv[1:] or ["A_v3", "B_zp", "C_zp_mix", "D_zp_mix_nearest", "E_zp_unpaired", "D_seed2"]
T = {v: table(v) for v in V}
print("smooth best age per run [Gyr]  (grid-minimum age in brackets)")
runs = T[V[0]].index
print("run    " + "".join(f"{v.split('_')[0] + ('2' if 'seed2' in v else ''):>16s}" for v in V))
for r in runs:
    print(f"{r}  " + "".join(f"{T[v].loc[r, 'smooth_age']:8.2f} ({T[v].loc[r, 'grid_age']:4.1f})" if r in T[v].index else f"{'':>16s}" for v in V))
D = T["D_zp_mix_nearest"]
for v in V:
    if v == "D_zp_mix_nearest":
        continue
    t = T[v].join(D, rsuffix="_D", how="inner")
    ds, dg = t.smooth_age - t.smooth_age_D, t.grid_age - t.grid_age_D
    print(f"{v:18s} - D over {len(t)} runs: smooth-age shift median {ds.median():+.2f} (range {ds.min():+.2f}..{ds.max():+.2f}, "
          f"{int((ds < 0).sum())} younger / {int((ds > 0).sum())} older);  grid-age shift median {dg.median():+.2f}, "
          f"rms {np.sqrt((dg ** 2).mean()):.2f} vs smooth rms {np.sqrt((ds ** 2).mean()):.2f}")
