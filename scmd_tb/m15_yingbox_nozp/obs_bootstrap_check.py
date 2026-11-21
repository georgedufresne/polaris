#!/usr/bin/env python3
"""How big is the 2D KS from resampling the observed M15 stars themselves (a non-parametric bootstrap),
compared with our null (model sCMD vs model sCMD at the best cell: 0.0076 +/- 0.0026)?

  obs_bootstrap_check.py SETUP_JSON N_RESAMPLES

Scored exactly like the fits (stats.ks_2d, 95th percentile, KS window, the observed fit stars of SETUP):
  (a) the observed stars against a resample of themselves with replacement (the observed set is the
      reference, the resample plays the sCMD);
  (b) two independent resamples against each other;
  (c) two random halves of the observed stars against each other (half the stars each).
"""
import json
import sys
import time

import numpy as np

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
sys.path.insert(0, f"{B}/m15_yingbox_nozp/code")
from isochronetoolbox.scmd import Setup, context                    # noqa: E402
from isochronetoolbox.scmd.stats import ks_2d                       # noqa: E402

d = {k: v for k, v in json.load(open(sys.argv[1])).items() if k in Setup.__dataclass_fields__}
n_res = int(sys.argv[2])
ctx = context(Setup(**d))
obs = ctx.obs3
N = len(obs)
rng = np.random.default_rng(20261111)
out = {"a": [], "b": [], "c": []}
t0 = time.time()
for _ in range(n_res):
    out["a"].append(ks_2d(obs, obs[rng.integers(0, N, N)], ctx.obs3_cdf)[0])
    r1, r2 = obs[rng.integers(0, N, N)], obs[rng.integers(0, N, N)]
    out["b"].append(ks_2d(r1, r2)[0])
    p = rng.permutation(N)
    out["c"].append(ks_2d(obs[p[:N // 2]], obs[p[N // 2:]])[0])
print(f"{N} observed stars in the KS window, {n_res} resamples each ({time.time() - t0:.0f} s)")
names = {"a": "observed vs a resample of itself", "b": "two independent resamples",
         "c": "two random halves (N/2 each)"}
for k, v in out.items():
    v = np.array(v)
    print(f"  ({k}) {names[k]:34s} KS mean {v.mean():.5f}  std {v.std():.5f}  mean + 3 sigma {v.mean() + 3 * v.std():.5f}")
