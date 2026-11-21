#!/usr/bin/env python3
"""Checks before the <cluster>_yingbox_nozp grid.

  validate_yingbox.py m15|tuc

1. With the production setup, this experiment's code copy reproduces the production best refined cell
   (<cluster>_final_refined) to the KS in the file.
2. With this experiment's setup (zero-point map off; for M15 also E 0.08-0.15), the same cell: KS, number of
   fit stars, and that no zero-point map is applied.
"""
import json
import sys
import time

import numpy as np
import pandas as pd

X = sys.argv[1]
B = "/dartfs-hpc/scratch/f007058/scmd_tb"
O = f"{B}/{X}_yingbox_nozp"
sys.path.insert(0, f"{O}/code")
from isochronetoolbox.scmd import Setup, context                    # noqa: E402

LIM = {"m15": "/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json",
       "tuc": f"{B}/tuc_final/tuc_age_limits_d4.json"}[X]
lim = {int(k): v for k, v in json.load(open(LIM)).items()}
g = pd.read_csv(f"{B}/{X}_final_refined/combined_2d_ks_results.csv")
g = g[g.age_gyr < g.simulation_id.map(lim) - 1e-6]
row = g.loc[g.ks_statistic.idxmin()]
run, age, mu, e = int(row.simulation_id), float(row.age_gyr), float(row.distance_modulus), float(row.extinction_ebv)


def ks_at(setup_path, override):
    d = {k: v for k, v in json.load(open(setup_path)).items() if k in Setup.__dataclass_fields__}
    d.update(override)
    S = Setup(**d)
    t0 = time.time()
    ctx = context(S)
    _, make = ctx.model(run, age)
    out = ctx.stage3(make(np.random.RandomState(S.crn_seed), mu, e))
    return S, ctx, out, time.time() - t0


S_p, ctx_p, out_p, t = ks_at(f"{B}/{X}_final/setup.json", json.load(open(f"{B}/{X}_final/refine/setup_override.json")))
print(f"{X} 1. production cell run {run}, {age} Gyr, mu {mu}, E {e}: file KS {row.ks_statistic:.6f}, "
      f"this code with the production setup {out_p['ks_statistic']:.6f} "
      f"({'OK' if abs(out_p['ks_statistic'] - row.ks_statistic) < 1e-9 else 'MISMATCH'}), fit stars {ctx_p.n_target} ({t:.0f} s)")
S_n, ctx_n, out_n, t = ks_at(f"{O}/setup.json", json.load(open(f"{O}/refine/setup_override.json")))
print(f"{X} 2. this setup: zp_correct {S_n.zp_correct}, zero-point map applied: {ctx_n.zp_map is not None}; "
      f"E grid {S_n.ebv_grid}; fit stars {ctx_n.n_target} (production {ctx_p.n_target}); "
      f"same cell KS {out_n['ks_statistic']:.6f} (production {out_p['ks_statistic']:.6f}) ({t:.0f} s)")
