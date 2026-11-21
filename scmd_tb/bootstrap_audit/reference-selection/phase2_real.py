"""Phase 2: realisations of the best cell scored against every reference (paired design).

usage: phase2_real.py split|shared J0 J1
  split : realisations j = J0..J1-1 with ctx.ast = half_b (np.random.RandomState(j), as the
          published split null) against R0..R9 (half_a references)
  shared: realisations with the full AST list against the data, S0..S5 and G0 (the grid's
          own selected sCMD); j as the published shared null
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/reference-selection")
from common import AGE, EVI, MU, OUT, RUN, Log, ks_many, load_ctx, prep_target, setup_path, window_points  # noqa: E402

setup_path("latest")
mode, j0, j1 = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
log = Log(f"phase2_{mode}_{j0}_{j1}")
S, ctx = load_ctx()
full = ctx.ast
if mode == "split":
    half_a, half_b = full.halves()
    ctx.ast = half_b
    names = [f"R{k}" for k in range(10)]
else:
    names = ["data"] + [f"S{k}" for k in range(6)] + ["G0"]
targets = []
for nm in names:
    if nm == "data":
        targets.append(prep_target(ctx.obs3, ctx.obs3_cdf))
    else:
        targets.append(prep_target(np.load(f"{OUT}/ref_{nm}.npy")))
log("targets", names, [len(t[0]) for t in targets])
_, make = ctx.model(RUN, AGE)
rows = []
for j in range(j0, j1):
    p = window_points(ctx, S, make(np.random.RandomState(j), MU, EVI))
    res = ks_many(p, targets)
    row = {"seed": j, "n": len(p)}
    for nm, r in zip(names, res):
        for k, v in r.items():
            row[f"{nm}_{k}"] = v
    rows.append(row)
    log(f"realisation {j}: " + " ".join(f"{nm} {r['ks']:.6f}" for nm, r in zip(names, res)))
pd.DataFrame(rows).to_csv(f"{OUT}/phase2_{mode}_{j0:03d}_{j1:03d}.csv", index=False)
log("done")
