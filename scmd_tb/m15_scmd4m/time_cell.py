#!/usr/bin/env python3
"""Time one cell (production best: run 10357, 14.2 Gyr, mu 15.58, E 0.096) with model sCMDs of the given sizes.
  time_cell.py SETUP_JSON N1 [N2 ...]   (N = 0: observed-sized, as production)"""
import json, sys, time, resource
import numpy as np
sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/m15_scmd4m/code")
from isochronetoolbox.scmd import Setup, context
d = {k: v for k, v in json.load(open(sys.argv[1])).items() if k in Setup.__dataclass_fields__}
d.update({"ebv_grid": [0.0, 0.1205, 0.001]})
ctx = context(Setup(**d))
_, make = ctx.model(10357, 14.2)
for n in map(int, sys.argv[2:]):
    t0 = time.time(); cat = make(np.random.RandomState(ctx.s.crn_seed), 15.58, 0.096, n=n or None); t1 = time.time()
    out = ctx.stage3(cat); t2 = time.time()
    print(f"n {len(cat):>8d}: generate {t1 - t0:6.1f} s, KS {t2 - t1:6.1f} s, KS {out['ks_statistic']:.6f} "
          f"(n_sim in window {out['n_simulated']}), max RSS {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6:.1f} GB", flush=True)
