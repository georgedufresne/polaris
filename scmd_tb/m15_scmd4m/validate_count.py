#!/usr/bin/env python3
"""Checks of stats.ecdf_at / ks_2d_count (the exact-count model ECDF for large sCMDs).
1. ecdf_at equals a brute-force count on random point sets (ties included).
2. On the production best cell (run 10357, 14.2 Gyr, mu 15.58, E 0.096), the KS with the interpolated and the
   counted model ECDF for model sCMDs of 32k, 250k and 1M stars, with timings."""
import json, sys, time
import numpy as np
sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/m15_scmd4m/code")
from isochronetoolbox.scmd import Setup, context
from isochronetoolbox.scmd.stats import ecdf_at, ks_2d, ks_2d_count
rng = np.random.default_rng(1)
for M, N, blk in ((5000, 300, 64), (20000, 1000, 4096), (3000, 500, 7)):
    sim, q = rng.random((M, 2)), rng.random((N, 2))
    sim[:50] = q[:50]                                    # exact ties
    bf = ((sim[None, :, 0] < q[:, None, 0]) & (sim[None, :, 1] < q[:, None, 1])).sum(1) / M
    print(f"1. M {M} N {N} block {blk}: max |count - brute force| {np.abs(ecdf_at(sim, q, blk) - bf).max():.2e}")
d = {k: v for k, v in json.load(open(sys.argv[1])).items() if k in Setup.__dataclass_fields__}
d.update({"ebv_grid": [0.0, 0.1205, 0.001]})
ctx = context(Setup(**d))
_, make = ctx.model(10357, 14.2)
for n in (0, 250000, 1000000):
    cat = make(np.random.RandomState(ctx.s.crn_seed), 15.58, 0.096, n=n or None)
    sim = ctx.points(cat)
    t0 = time.time(); k1, d1 = ks_2d(ctx.obs3, sim, ctx.obs3_cdf, 95.0); t1 = time.time()
    k2, d2 = ks_2d_count(ctx.obs3, sim, ctx.obs3_cdf, 95.0); t2 = time.time()
    print(f"2. n {len(sim):>8d}: KS interp {k1:.6f} ({t1 - t0:5.1f} s)  count {k2:.6f} ({t2 - t1:5.1f} s)  "
          f"per-star |interp - count| median {np.median(np.abs(d1 - d2)):.1e}, 99% {np.percentile(np.abs(d1 - d2), 99):.1e}", flush=True)
