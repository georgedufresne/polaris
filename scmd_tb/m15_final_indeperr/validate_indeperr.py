#!/usr/bin/env python3
"""Checks before the m15_final_indeperr grid (production with independent per-filter artificial-star errors).

1. This code copy with ../m15_final_refined/setup.json reproduces production's best cell (run 10357, 14.2 Gyr,
   mu 15.58, E 0.096) to the KS in its file.
2. The same cell with ast_errors "independent" (this experiment's setup + the refinement E grid): KS, fit-star and
   sCMD counts (crn seed, one sCMD each). Also run 10392 at 13.3 Gyr, mu 15.51, E 0.065 (the error-chain cell)."""
import json
import sys
import time

import numpy as np
import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
H = f"{B}/m15_final_indeperr"
sys.path.insert(0, f"{H}/code")
from isochronetoolbox.scmd import Setup, context    # noqa: E402


def setup(base, override):
    d = {k: v for k, v in json.load(open(base)).items() if k in Setup.__dataclass_fields__}
    d.update(override)
    return context(Setup(**d))


def ks(ctx, run, age, mu, e):
    _, make = ctx.model(run, age)
    o = ctx.stage3(make(np.random.RandomState(ctx.s.crn_seed), mu, e))
    return o["ks_statistic"], o["n_observed"], o["n_simulated"]


t0 = time.time()
ref = pd.read_csv(f"{B}/m15_final_refined/combined_2d_ks_results.csv")
r = ref[(ref.simulation_id == 10357) & np.isclose(ref.age_gyr, 14.2) & np.isclose(ref.distance_modulus, 15.58)
        & np.isclose(ref.extinction_ebv, 0.096)].ks_statistic.min()
P = f"{B}/m15_final_refined/setup.json"
prod = setup(P, {})
k, no, ns = ks(prod, 10357, 14.2, 15.58, 0.096)
print(f"1. production setup, this code: KS {k:.7f} vs file {r:.7f} (diff {k - r:+.1e}); {no} fit, {ns} sCMD "
      f"[{time.time() - t0:.0f} s]", flush=True)
over = json.load(open(f"{H}/refine/setup_override.json"))
ind = setup(f"{H}/setup.json", over)
print(f"   this experiment: ast_errors {ind.s.ast_errors}, radial {ind.s.radial}, martin_sign {ind.s.martin_sign}, "
      f"ks_percentile {ind.s.ks_percentile}, zp_correct {ind.s.zp_correct}")
for cell in [(10357, 14.2, 15.58, 0.096), (10392, 13.3, 15.51, 0.065)]:
    a = ks(prod, *cell)
    b = ks(ind, *cell)
    print(f"2. {cell}: KS joint (production) {a[0]:.5f}, independent {b[0]:.5f}; sCMD stars {a[2]} / {b[2]} "
          f"[{time.time() - t0:.0f} s]", flush=True)
