#!/usr/bin/env python3
"""Colour width about the observed ridge, observed vs model, by radius: before and after the zero-point map."""
import json, sys
import numpy as np
E = "/dartfs-hpc/scratch/f007058/scmd_tb/m15_errchain"
sys.path.insert(0, f"{E}/code")
from isochronetoolbox.scmd import Setup, context
base = {k: v for k, v in json.load(open(f"{E}/setup.json")).items() if k in Setup.__dataclass_fields__}
base.update(imf_slope=-1.0, binary_fraction=0.017, ast_errors="nearest", radial="corrected", radial_source="raw")
RB = np.array([0, 100, 200, 400, 700, 1200, 2000, 4500])
def rsig(x):
    return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) > 50 else np.nan
for name, kw, cell in (("no zero-point map (V3)", {}, (10392, 13.3, 15.51, 0.065)),
                       ("zero-point map, both bands (V6)", dict(zp_correct=True, zp_ref=(19.3, 20.6), zp_mag_frac=0.57),
                        (10392, 13.6, 15.475, 0.06))):
    d = dict(base); d.update(kw); S = Setup(**d); C = context(S)
    fit = C.ast._args[0]; xc, yc = C.ast.centre_xy
    fr = np.hypot(fit.x.values - xc, fit.y.values - yc); fv, fc = fit.v.values, (fit.v - fit.i).values
    cat = C.model(cell[0], cell[1])[1](np.random.RandomState(S.crn_seed), cell[2], cell[3], with_radius=True)
    rid = C.rcut.ridge; m1 = rid.mag > 19.3; o1 = np.argsort(rid.mag[m1])
    f = lambda V: np.interp(V, rid.mag[m1][o1], rid.color[m1][o1])
    print(f"\n{name}, model run {cell[0]} {cell[1]} Gyr mu {cell[2]} E {cell[3]}")
    print("  MS 19.4-20.6   r (px)     obs     model   ratio   missing sigma (quadrature)")
    for k in range(len(RB) - 1):
        so = (fv > 19.4) & (fv < 20.6) & (fr >= RB[k]) & (fr < RB[k + 1])
        sm = (cat[:, 0] > 19.4) & (cat[:, 0] < 20.6) & (cat[:, 2] >= RB[k]) & (cat[:, 2] < RB[k + 1])
        a = rsig(fc[so] - f(fv[so])); b = rsig((cat[sm, 0] - cat[sm, 1]) - f(cat[sm, 0]))
        print(f"            {RB[k]:5d}-{RB[k+1]:5d}  {a:.4f}  {b:.4f}   {a / b:.2f}    {np.sqrt(max(a*a - b*b, 0)):.4f}")
