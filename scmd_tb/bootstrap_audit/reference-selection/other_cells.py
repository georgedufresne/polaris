"""Data vs fresh realisations (full AST list, latest code) of the grid's next-best cells:
their expected KS and per-cell scatter, against their grid values.

usage: other_cells.py N
"""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/reference-selection")
from common import OUT, Log, ks_many, load_ctx, prep_target, setup_path, window_points  # noqa: E402

setup_path("latest")
N = int(sys.argv[1])
CELLS = [(10005, 11.1, 13.43, 0.030), (10091, 11.4, 13.42, 0.035),
         (10099, 10.9, 13.42, 0.020), (10030, 12.3, 13.34, 0.005)]
log = Log("other_cells")
S, ctx = load_ctx()
grid = pd.read_csv("/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02/combined_2d_ks_results.csv")
data = prep_target(ctx.obs3, ctx.obs3_cdf)
rows = []
for run, age, mu, e in CELLS:
    g = grid[(grid.simulation_id == run) & np.isclose(grid.age_gyr, age) & np.isclose(grid.distance_modulus, mu)
             & np.isclose(grid.extinction_ebv, e)].ks_statistic.values[0]
    _, make = ctx.model(run, age)
    ks = []
    for j in range(N):
        p = window_points(ctx, S, make(np.random.RandomState(j), mu, e))
        r = ks_many(p, [data])[0]
        rows.append({"run": run, "age": age, "mu": mu, "evi": e, "grid_ks": g, "seed": j, "n": len(p),
                     **{f"data_{k}": v for k, v in r.items()}})
        ks.append(r["ks"])
    log(f"cell {run} {age} {mu} {e}: grid {g:.6f}; fresh mean {np.mean(ks):.6f} std {np.std(ks, ddof=1):.6f} "
        f"(n {N}); grid - mean = {g - np.mean(ks):+.6f}")
    pd.DataFrame(rows).to_csv(f"{OUT}/other_cells.csv", index=False)
log("done")
