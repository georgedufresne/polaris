#!/usr/bin/env python3
"""Checks before the m15_yingiso_nozp grid (Ying's published isochrones; otherwise m15_yingbox_nozp: no zero-point map,
E 0.08-0.15 applied directly, our sCMDs and our 95th-percentile KS).

1. This experiment's code copy with m15_yingbox_nozp's refined setup (our isochrones) reproduces that run's best cell
   (run 10245, 14.0 Gyr, mu 15.54, E 0.085) to the KS in its file.
2. The same cell with Ying's isochrones (this setup): KS, fit-star and sCMD counts.
3. Ying's published M15 best fit (DM 15.49, E 0.08) for runs 10245, 10357, 10392 at 13.4 and 14.2 Gyr: KS with his
   isochrones and with ours (crn seed, one sCMD each)."""
import json
import sys
import time

import numpy as np
import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
H = f"{B}/m15_yingiso_nozp"
sys.path.insert(0, f"{H}/code")
from isochronetoolbox.scmd import Setup, context    # noqa: E402


def ctx_for(path):
    d = {k: v for k, v in json.load(open(path)).items() if k in Setup.__dataclass_fields__}
    return context(Setup(**d))


def ks(ctx, run, age, mu, e):
    _, make = ctx.model(run, age)
    o = ctx.stage3(make(np.random.RandomState(ctx.s.crn_seed), mu, e))
    return o["ks_statistic"], o["n_observed"], o["n_simulated"]


t0 = time.time()
ref = pd.read_csv(f"{B}/m15_yingbox_nozp_refined/combined_2d_ks_results.csv")
r = ref[(ref.simulation_id == 10245) & np.isclose(ref.age_gyr, 14.0) & np.isclose(ref.distance_modulus, 15.54)
        & np.isclose(ref.extinction_ebv, 0.085)].iloc[0]
ours = ctx_for(f"{B}/m15_yingbox_nozp_refined/setup.json")
k, no, ns = ks(ours, 10245, 14.0, 15.54, 0.085)
print(f"1. our isochrones, this code copy: KS {k:.7f} vs file {r.ks_statistic:.7f} (diff {k - r.ks_statistic:+.1e}); "
      f"{no} fit stars, {ns} sCMD stars  [{time.time() - t0:.0f} s]", flush=True)
his = ctx_for(f"{H}/setup.json")
k2, no2, ns2 = ks(his, 10245, 14.0, 15.54, 0.085)
print(f"2. Ying's isochrones, same cell: KS {k2:.5f}; {no2} fit stars, {ns2} sCMD stars; zp_correct {his.s.zp_correct}, "
      f"iso_format {his.s.iso_format}  [{time.time() - t0:.0f} s]", flush=True)
print("3. at Ying's published best fit (DM 15.49, E 0.08): KS with his isochrones / ours")
for run in (10245, 10357, 10392):
    for age in (13.4, 14.2):
        print(f"   run {run} {age} Gyr: {ks(his, run, age, 15.49, 0.08)[0]:.4f} / {ks(ours, run, age, 15.49, 0.08)[0]:.4f}",
              flush=True)
print(f"({time.time() - t0:.0f} s)")
