#!/usr/bin/env python3
"""Does a bigger model sCMD (2x, 4x, 16x the observed count) change the KS, its noise, and the age
discrimination? And what does duplicating the observed stars do?

For the cells that decide the age (production: 14.2 vs 12.3 Gyr; Ying's box without the map: 14.0 vs 12.3 Gyr),
the KS with model sCMDs of 1x, 2x, 4x, 16x the observed count, each with several seeds (crn_seed + k), the
model ECDF counted exactly (stats.ecdf_at) at every size so the sizes are comparable. Output in
logs/model_size_check.out."""
import json
import sys
import time

import numpy as np

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
sys.path.insert(0, f"{B}/m15_scmd4m/code")
from isochronetoolbox.scmd import Setup, context                    # noqa: E402
from isochronetoolbox.scmd.stats import ks_2d_count, rank2d         # noqa: E402

CELLS = {"m15_final_refined": [(10357, 14.2, 15.58, 0.096, "best, 14.2 Gyr"), (10392, 12.3, 15.60, 0.080, "12.3 Gyr")],
         "m15_yingbox_nozp_refined": [(10245, 14.0, 15.54, 0.085, "best, 14.0 Gyr"), (10392, 12.3, 15.60, 0.080, "12.3 Gyr")]}
SIZES = [(1, 5), (2, 5), (4, 5), (16, 3)]          # (multiple of the observed count, seeds)
t00 = time.time()
for base, cells in CELLS.items():
    d = {k: v for k, v in json.load(open(f"{B}/{base}/setup.json")).items() if k in Setup.__dataclass_fields__}
    d["ks_model_ecdf"] = "count"
    ctx = context(Setup(**d))
    N = ctx.n_target
    print(f"\n{base} ({N} observed fit stars)")
    res = {}
    for run, age, mu, e, name in cells:
        _, make = ctx.model(run, age)
        for f, ns in SIZES:
            ks = [ctx.stage3(make(np.random.RandomState(ctx.s.crn_seed + k), mu, e, n=f * N))["ks_statistic"]
                  for k in range(ns)]
            res[(name, f)] = np.array(ks)
            print(f"  {name:15s} model {f:2d}x ({f * N:>7d}): KS mean {np.mean(ks):.5f}  seed-to-seed std {np.std(ks):.5f}  "
                  f"(seed 0: {ks[0]:.5f})", flush=True)
    names = [c[4] for c in cells]
    print("  age discrimination, KS(12.3) - KS(best), mean over seeds (seed std of the difference):")
    for f, ns in SIZES:
        diff = res[(names[1], f)] - res[(names[0], f)]
        print(f"    model {f:2d}x: {diff.mean():+.5f} ({diff.std():.5f})")
    if base == "m15_final_refined":
        # duplicating the observed stars: every star counted twice
        run, age, mu, e, _ = cells[0]
        _, make = ctx.model(run, age)
        sim = ctx.points(make(np.random.RandomState(ctx.s.crn_seed), mu, e, n=N))
        o2 = np.vstack([ctx.obs3, ctx.obs3])
        k1 = ks_2d_count(ctx.obs3, sim, ctx.obs3_cdf)[0]
        k2 = ks_2d_count(o2, sim, rank2d(o2) / len(o2))[0]
        print(f"  observed stars duplicated (2 x {len(ctx.obs3)}): KS {k2:.6f} vs {k1:.6f} with each star once")
print(f"\n({time.time() - t00:.0f} s)")
