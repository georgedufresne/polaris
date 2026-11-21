#!/usr/bin/env python3
"""What age does the production pipeline give if the distance is fixed at Gaia's?

The grid's mu is the apparent F606W distance modulus (reddening "evi": F606W = M + mu, F814W = M + mu - E), so a
true modulus mu0 fixes mu = mu0 + A_F606W = mu0 + 2.82 E (ACS: A_F606W 2.88 E(B-V), E(F606W-F814W) 1.02 E(B-V)).
Baumgardt & Vasiliev (2021): 10.71 +/- 0.10 kpc, mu0 = 15.149 +/- 0.020.

For a few runs and ages, scan E along that line (mu0 = 15.149, and mu0 = 15.189 = Gaia + 2 sigma, toward where
the free fits sit) with the production setup unchanged (zero-point map on, 32k sCMDs, crn seed, stage-3 KS). The
model ECDF is counted exactly (stats.ecdf_at), which differs from the interpolated production value by about
0.0002. Output: logs/gaia_line_scan.out and gaia_line_scan.csv."""
import json
import sys
import time

import numpy as np
import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
H = f"{B}/m15_gaia_distance"
sys.path.insert(0, f"{B}/m15_scmd4m/code")       # production code + the counted ECDF option (n_model 0 = production)
from isochronetoolbox.scmd import Setup, context  # noqa: E402

A606 = 2.88 / 1.02
MU0 = [15.149, 15.189]
RUNS = [10357, 10392, 10245]                      # production best; old 13.3 Gyr escan run; no-map best
AGES = [12.6, 13.0, 13.4, 13.8, 14.2, 14.6]
EGRID = np.round(np.arange(0.06, 0.1601, 0.01), 3)
lim = {int(k): v for k, v in json.load(open("/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json")).items()}

d = {k: v for k, v in json.load(open(f"{B}/m15_final_refined/setup.json")).items() if k in Setup.__dataclass_fields__}
d["ks_model_ecdf"] = "count"
ctx = context(Setup(**d))
free = pd.read_csv(f"{B}/m15_final_refined/combined_2d_ks_results.csv",
                   usecols=["simulation_id", "age_gyr", "distance_modulus", "extinction_ebv", "ks_statistic"])
t0 = time.time()
rows = []
for run in RUNS:
    for age in AGES:
        if age >= lim[run] - 1e-6:
            continue
        f = free[(free.simulation_id == run) & np.isclose(free.age_gyr, age)]
        fb = f.loc[f.ks_statistic.idxmin()]
        _, make = ctx.model(run, age)
        line = []
        for mu0 in MU0:
            for e in EGRID:
                mu = mu0 + A606 * e
                ks = ctx.stage3(make(np.random.RandomState(ctx.s.crn_seed), mu, e))["ks_statistic"]
                rows.append(dict(run=run, age=age, mu0=mu0, E=e, mu=round(mu, 4), ks=ks))
                line.append((mu0, e, ks))
        best = {m: min((x for x in line if x[0] == m), key=lambda x: x[2]) for m in MU0}
        print(f"run {run} {age:4.1f} Gyr | free fit: KS {fb.ks_statistic:.4f} at mu {fb.distance_modulus:.2f} E {fb.extinction_ebv:.3f}"
              f" (mu0 {fb.distance_modulus - A606 * fb.extinction_ebv:.3f}) | "
              + " | ".join(f"Gaia mu0 {m}: KS {best[m][2]:.4f} at E {best[m][1]:.2f}" for m in MU0)
              + f"  [{time.time() - t0:.0f} s]", flush=True)
pd.DataFrame(rows).to_csv(f"{H}/gaia_line_scan.csv", index=False)
r = pd.DataFrame(rows)
print("\nbest KS along the Gaia line, by age (minimum over the runs above):")
print(r.groupby(["mu0", "age"]).ks.min().unstack(0).round(4).to_string())
