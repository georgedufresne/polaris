#!/usr/bin/env python3
"""Check the KS descent against exhaustive scoring: the error-model chains (m15_errchain, tuc_errchain,
variant V6_zpsplit = the final setups) scored every cell of a local (age, mu, E) grid for 8 runs. For every
(run, age) of those runs that the refinement has finished, the refined minimum must be at or below the
chain's minimum over the cells the two grids share (same catalogues: crn, identical setups).

  validate.py m15|tuc
"""
import glob
import os
import sys

import numpy as np
import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
C = {"m15": ("m15_final", "m15_errchain"), "tuc": ("tuc_final", "tuc_errchain")}[sys.argv[1]]
G, E = f"{B}/{C[0]}", f"{B}/{C[1]}"
ch = pd.concat([pd.read_csv(p) for p in glob.glob(f"{E}/fits/V6_zpsplit/*.csv")], ignore_index=True)
grid = pd.concat([pd.read_csv(p) for r in ch.run.unique() for p in glob.glob(f"{G}/units/{r}/*.ks.csv")],
                 ignore_index=True)
ref = [pd.read_csv(p) for r in ch.run.unique() for p in glob.glob(f"{G}/refine/units/{r}/*.csv")]
ref = pd.concat([x for x in ref if len(x)], ignore_index=True) if ref else grid.iloc[:0]
done_runs = {int(os.path.basename(os.path.dirname(p))) for p in glob.glob(f"{G}/refine/units/*/*.csv")}
m = pd.concat([grid, ref], ignore_index=True)
dm = np.round(np.sort(m.distance_modulus.unique()), 3)
ev = np.round(np.sort(m.extinction_ebv.unique()), 3)
shared = ch[ch.mu.round(3).isin(dm) & ch.evi.round(3).isin(ev)]
cm = shared.groupby(["run", "age"]).ks.min()
gm = grid.groupby(["simulation_id", "age_gyr"]).ks_statistic.min()
mm = m.groupby(["simulation_id", "age_gyr"]).ks_statistic.min()
rows = []
for (r, a), v in cm.items():
    if (r, a) not in mm.index or r not in done_runs:
        continue
    rows.append((r, a, v, gm.get((r, a), np.nan), mm[(r, a)]))
t = pd.DataFrame(rows, columns=["run", "age", "chain_shared_min", "grid_min", "refined_min"])
if t.empty:
    sys.exit("no refined (run, age) of the chain runs yet")
bad = t[t.refined_min > t.chain_shared_min + 1e-9]
print(f"{C[0]}: {len(t)} (run, age) of the chain runs refined; stage 2 + 3 alone missed the shared-cell minimum "
      f"in {(t.grid_min > t.chain_shared_min + 1e-9).sum()}; refinement misses it in {len(bad)}")
print(f"  median KS: grid {t.grid_min.median():.5f}, refined {t.refined_min.median():.5f}, chain shared cells "
      f"{t.chain_shared_min.median():.5f}")
if len(bad):
    print(bad.assign(excess=bad.refined_min - bad.chain_shared_min).sort_values("excess", ascending=False)
          .head(15).to_string(index=False, float_format=lambda x: f"{x:.5f}"))
