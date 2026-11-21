#!/usr/bin/env python3
"""Checks before the m15_final_yingiso_indeperr grid (production with Ying's isochrones and independent per-filter errors).

1. This code with ../m15_yingiso_nozp_refined/setup.json reproduces that run's best cell (10089, 12.6 Gyr, mu 15.50,
   E 0.086) to the KS in its file. The file's E comes from summed refinement steps, so a rounding-level difference is
   expected, as in ../m15_yingiso_err99.
2. KS at a few cells (one sCMD each, crn seed) for two setups: production with Ying's isochrones and production's joint
   errors (one change), and this setup (+ independent errors). Cells: that best cell, production's best cell (10357,
   14.2, 15.58, 0.096) and the replica's best cells (10411, 13.2, 15.429, 0.080; 10336, 13.2, 15.443, 0.085).
Output: logs/validate_yingiso_indeperr.out."""
import json
import sys
import time

import numpy as np
import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
H = f"{B}/m15_final_yingiso_indeperr"
sys.path.insert(0, f"{H}/code")
from isochronetoolbox.scmd import Setup, context    # noqa: E402


def ctx_of(base, override):
    d = {k: v for k, v in json.load(open(base)).items() if k in Setup.__dataclass_fields__}
    d.update(override)
    return context(Setup(**d))


def ks(ctx, run, age, mu, e):
    _, make = ctx.model(run, age)
    o = ctx.stage3(make(np.random.RandomState(ctx.s.crn_seed), mu, e))
    return o["ks_statistic"]


t0 = time.time()
ref = pd.read_csv(f"{B}/m15_yingiso_nozp_refined/combined_2d_ks_results.csv")
r = ref[(ref.simulation_id == 10089) & np.isclose(ref.age_gyr, 12.6) & np.isclose(ref.distance_modulus, 15.50)
        & np.isclose(ref.extinction_ebv, 0.086)].ks_statistic.min()
k = ks(ctx_of(f"{B}/m15_yingiso_nozp_refined/setup.json", {}), 10089, 12.6, 15.50, 0.086)
print(f"1. ../m15_yingiso_nozp setup, this code: KS {k:.7f} vs file {r:.7f} (diff {k - r:+.1e}) [{time.time() - t0:.0f} s]",
      flush=True)
over = json.load(open(f"{H}/refine/setup_override.json"))
joint = ctx_of(f"{H}/setup.json", {**over, "ast_errors": "nearest"})
ind = ctx_of(f"{H}/setup.json", over)
print(f"   this setup: iso {ind.s.iso.split('/')[-1]}, iso_format {ind.s.iso_format}, ast_errors {ind.s.ast_errors}, "
      f"zp_correct {ind.s.zp_correct}, ks_percentile {ind.s.ks_percentile}, E grid {ind.s.ebv_grid}")
print("2. KS: production + Ying's isochrones (joint errors) | + independent errors (this setup)")
for cell in [(10089, 12.6, 15.50, 0.086), (10357, 14.2, 15.58, 0.096), (10411, 13.2, 15.429, 0.080),
             (10336, 13.2, 15.443, 0.085)]:
    print(f"   {cell}: {ks(joint, *cell):.5f} | {ks(ind, *cell):.5f}  [{time.time() - t0:.0f} s]", flush=True)
