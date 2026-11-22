#!/usr/bin/env python3
"""Is the W v1.4 NGC 6791 age bias the [Fe/H] -> Z convention?

Claim to test (from a separate analysis of W's training grid): W's training tracks have Z = 1.121 x getz([Fe/H], [a/Fe]),
whereas our DSEP3 (and OG) runs use Z = getz([Fe/H], [a/Fe]) with the same gs98getz (checked: the isochrone headers of
OG, DSEP3 and the var files agree). gen_W14_iso.py passes the var file's [Fe/H] straight to W, so W would model 12 % more
metals than the DSEP3 run of the same var file.

For each run, two W v1.4 central (mode) isochrone grids, 7.0-9.9 Gyr:
  nominal:   [Fe/H] as written in the var file (what the W14 sCMD leg used);
  converted: [Fe/H]_W = getz^-1( getz([Fe/H]) / 1.121 ) at the run's [a/Fe] (age head and flow both get it).
The run's DSEP3 physics-on isochrone at 8.0 / 8.3 / 8.6 Gyr is fitted against each grid (same fitter as
../flow_noise: turnoff box, IMF weights, free distance + reddening shift). If the claim holds, "converted" should give
W age - DSEP3 age near 0 and remove the ~0.002 dex turnoff temperature offset; "nominal" should reproduce the bias.

  python feh_test.py RUN [RUN ...] [--factor 1.121]
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np
import torch

S = "/dartfs-hpc/scratch/f007058/emul_W14"
for p in (f"{S}/flow_noise", S, f"{S}/CONF1DENCE-W_fcfac3c"):
    sys.path.insert(0, p)
from flow_noise import GRID, MASS_HI, MASS_LO, NMASS, Model, dsep3_block, fit_age, pseudo, turnoff   # noqa: E402
from gen_W14_iso import age_table, invert, phys_for                                               # noqa: E402
from isochrone_W import central_mode                                                              # noqa: E402
from predict_W import AgeHead, WEmulator                                                          # noqa: E402


def getz(feh, afe):
    return float(subprocess.run([f"{S}/gs98getz", f"{feh:.6f}", f"{afe:.6f}"], capture_output=True, text=True).stdout)


def feh_for_z(z, afe, lo=-2.0, hi=1.5):
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if getz(mid, afe) < z:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+", type=int)
    ap.add_argument("--factor", type=float, default=1.121)
    ap.add_argument("--out", default=f"{S}/feh_convention/output")
    a = ap.parse_args()
    torch.set_num_threads(int(os.environ.get("SLURM_CPUS_PER_TASK", "4")))
    os.makedirs(a.out, exist_ok=True)
    W = WEmulator(); ah = AgeHead()
    lh = W.fmt["label_header"]; default = W.fmt["label_default"]
    for run in a.runs:
        out = f"{a.out}/run_{run}.json"
        if os.path.exists(out):
            print(f"{run}: exists", flush=True); continue
        t0 = time.time()
        phys = phys_for(run, f"{S}/var", lh, default, "var")
        feh, afe = phys[lh.index("FeH")], phys[lh.index("AlphaFe")]
        zr = getz(feh, afe); fehc = feh_for_z(zr / a.factor, afe)
        res = {"run": run, "factor": a.factor, "feh": feh, "afe": afe, "z_run": zr, "feh_converted": fehc,
               "z_check": getz(fehc, afe) * a.factor}
        for label, f in (("nominal", feh), ("converted", fehc)):
            p = list(phys); p[lh.index("FeH")] = f
            p22 = np.asarray(p[:22], float)[:len(ah.theta)]
            eg, mg, LA = age_table(ah, p22, MASS_LO, MASS_HI, NMASS)
            cen, models = {}, []
            for age in GRID:
                m, e = invert(ah, p22, eg, mg, LA, float(age), MASS_LO, MASS_HI)
                o = central_mode(W, e, m, p)
                x, y = pseudo(o["log_Teff"], o["log_L"])
                cen[float(age)] = (x, y, o["log_L"])
                models.append(Model(float(age), x, y, o["log_L"]))
            res[label] = {}
            for age in (8.0, 8.3, 8.6):
                lT, lL, mm = dsep3_block(run, age); x, y = pseudo(lT, lL)
                ag, rms, ok = fit_age(models, x, y, lL, mm)
                xd, yd = turnoff(x, y, lL); xw, yw, lLw = cen[age]; xtw, ytw = turnoff(xw, yw, lLw)
                res[label][str(age)] = {"w_age_for_dsep3": ag, "rms": rms, "inside_grid": ok,
                                        "to_dlogTeff_d3_minus_w": float(-(xd - xtw) / 5.0),
                                        "to_dlogL_d3_minus_w": float(-(yd - ytw) / 2.5)}
        res["seconds"] = round(time.time() - t0, 1)
        json.dump(res, open(out + ".tmp", "w"), indent=1); os.replace(out + ".tmp", out)
        n = res["nominal"]["8.3"]; c = res["converted"]["8.3"]
        print(f"{run}: [Fe/H] {feh:+.3f} -> {fehc:+.3f} (Z {zr:.5f}, check {res['z_check']:.5f}); DSEP3 8.3 Gyr fits as "
              f"W nominal {n['w_age_for_dsep3']:.3f} / converted {c['w_age_for_dsep3']:.3f}; TO dlogTeff (D3-W) "
              f"{n['to_dlogTeff_d3_minus_w']:+.5f} / {c['to_dlogTeff_d3_minus_w']:+.5f}  ({res['seconds']:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
