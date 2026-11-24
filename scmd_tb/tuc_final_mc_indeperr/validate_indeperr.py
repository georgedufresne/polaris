#!/usr/bin/env python3
"""Checks before the tuc_final_mc_indeperr grid (the Monte Carlo-on 47 Tuc fit with independent per-filter
artificial-star errors, the advisor's / M. Ying's MakeTestCMD model). Same checks as ../m15_final_indeperr.

1. This code copy with ../tuc_final_mc/setup.json reproduces that grid's lowest-KS row to the KS in its file.
2. The same cell and two first-pass best cells (10567 10.2 Gyr mu 13.50 E 0.040; 10105 12.3 Gyr mu 13.35 E 0.010)
   with ast_errors "independent": KS, fit-star and sCMD counts (crn seed, one sCMD each)."""
import json
import sys
import time

import numpy as np
import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
H = f"{B}/tuc_final_mc_indeperr"
sys.path.insert(0, f"{H}/code")
from isochronetoolbox.scmd import Setup, context    # noqa: E402


def setup(base):
    d = {k: v for k, v in json.load(open(base)).items() if k in Setup.__dataclass_fields__}
    return context(Setup(**d))


def ks(ctx, run, age, mu, e):
    _, make = ctx.model(run, age)
    o = ctx.stage3(make(np.random.RandomState(ctx.s.crn_seed), mu, e))
    return o["ks_statistic"], o["n_observed"], o["n_simulated"]


t0 = time.time()
g = pd.read_csv(f"{B}/tuc_final_mc/combined_2d_ks_results.csv")
r = g.loc[g.ks_statistic.idxmin()]
cell = (int(r.simulation_id), float(r.age_gyr), float(r.distance_modulus), float(r.extinction_ebv))
prod = setup(f"{B}/tuc_final_mc/setup.json")
k, no, ns = ks(prod, *cell)
print(f"1. grid best {cell}: this code, joint errors KS {k:.7f} vs file {r.ks_statistic:.7f} (diff {k - r.ks_statistic:+.1e}); "
      f"{no} fit, {ns} sCMD [{time.time() - t0:.0f} s]", flush=True)
ind = setup(f"{H}/setup.json")
print(f"   this experiment: ast_errors {ind.s.ast_errors}, radial {ind.s.radial}, martin_sign {ind.s.martin_sign}, "
      f"zp_correct {ind.s.zp_correct}, cut {ind.s.cut}")
for c in [cell, (10567, 10.2, 13.50, 0.040), (10105, 12.3, 13.35, 0.010)]:
    a, b = ks(prod, *c), ks(ind, *c)
    print(f"2. {c}: KS joint {a[0]:.5f}, independent {b[0]:.5f}; sCMD stars {a[2]} / {b[2]} [{time.time() - t0:.0f} s]", flush=True)
