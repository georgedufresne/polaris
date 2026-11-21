#!/usr/bin/env python3
"""Checks before the m15_yingiso_err99 grid (Ying's isochrones + his artificial-star error model + 99th-percentile KS;
otherwise ../m15_yingiso_nozp).

1. This code copy with ../m15_yingiso_nozp_refined/setup.json reproduces that run's best cell (run 10089, 12.6 Gyr,
   mu 15.50, E 0.086) to the KS in its file.
2. At the same cell, the changes one at a time: + his error model (independent per-band draws from 10 x 12
   (radius, mag) bins, unweighted observed radii), + his error sign (observed = true + (input - output)), + the
   99th-percentile KS (= this setup). KS, fit-star and sCMD counts (crn seed, one sCMD each)."""
import json
import sys
import time

import numpy as np
import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
H = f"{B}/m15_yingiso_err99"
sys.path.insert(0, f"{H}/code")
from isochronetoolbox.scmd import Setup, context    # noqa: E402

RUN, AGE, MU, E = 10089, 12.6, 15.50, 0.086


def ks(base, override):
    d = {k: v for k, v in json.load(open(base)).items() if k in Setup.__dataclass_fields__}
    d.update(override)
    ctx = context(Setup(**d))
    _, make = ctx.model(RUN, AGE)
    o = ctx.stage3(make(np.random.RandomState(ctx.s.crn_seed), MU, E))
    return o["ks_statistic"], o["n_observed"], o["n_simulated"], ctx


t0 = time.time()
ref = pd.read_csv(f"{B}/m15_yingiso_nozp_refined/combined_2d_ks_results.csv")
r = ref[(ref.simulation_id == RUN) & np.isclose(ref.age_gyr, AGE) & np.isclose(ref.distance_modulus, MU)
        & np.isclose(ref.extinction_ebv, E)].ks_statistic.min()
P = f"{B}/m15_yingiso_nozp_refined/setup.json"
k, no, ns, _ = ks(P, {})
print(f"1. ../m15_yingiso_nozp setup, this code: KS {k:.7f} vs file {r:.7f} (diff {k - r:+.1e}); {no} fit, {ns} sCMD  "
      f"[{time.time() - t0:.0f} s]", flush=True)
steps = [("+ his error model (independent, observed radii)", {"ast_errors": "independent", "radial": "observed"}),
         ("+ his error sign", {"ast_errors": "independent", "radial": "observed", "martin_sign": True}),
         ("+ 99th-percentile KS", {"ast_errors": "independent", "radial": "observed", "martin_sign": True, "ks_percentile": 99.0})]
for name, ov in steps:
    k, no, ns, ctx = ks(P, ov)
    print(f"2. {name:48s} KS {k:.5f}; {no} fit, {ns} sCMD  [{time.time() - t0:.0f} s]", flush=True)
k, no, ns, ctx = ks(f"{H}/setup.json", {})
print(f"   this experiment's setup.json: KS {k:.5f}; ast_errors {ctx.s.ast_errors}, radial {ctx.s.radial}, martin_sign "
      f"{ctx.s.martin_sign}, ks_percentile {ctx.s.ks_percentile}, zp_correct {ctx.s.zp_correct}, iso_format {ctx.s.iso_format}")
