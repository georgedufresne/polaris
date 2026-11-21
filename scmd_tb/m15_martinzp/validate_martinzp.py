#!/usr/bin/env python3
"""Checks before the m15_martinzp grid (run from m15_martinzp/code).

  validate_martinzp.py M2_CATALOGUE M2_FIDUCIAL

1. ast.martin_zero_point reproduces M. Ying's own output: applied to the catalogue his
   M15_fitstars_ZPCR.dat was made from (M15_2023.dat; every one of its stars matches it, not our M15.dat),
   with V_SGB 18.85 as in M15_ZPCR.ipynb, it must give his corrected colours. Also against the notebook's
   brute-force neighbour loop on 300 stars of our M15.dat.
2. With zp_method "ours" (the default) this code copy reproduces production: the lowest-KS refined cell of
   m15_final_refined must come out with the same KS.
3. Our zero-point map against Martin's on our M15 data, and the KS of the same cell with Martin's map.
"""
import json
import sys
import time

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.spatial import cKDTree

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
sys.path.insert(0, f"{B}/m15_martinzp/code")
from isochronetoolbox.scmd import Setup, context                    # noqa: E402
from isochronetoolbox.scmd.ast import martin_zero_point              # noqa: E402

M2, M2FID = sys.argv[1], sys.argv[2]
MARTIN = "/dartfs-hpc/scratch/f007058/Files_from_Martin/M15_data"


def rawframe(path):
    r = pd.read_csv(path, sep=r"\s+")
    return pd.DataFrame({"v": r["Vvega"], "i": r["Ivega"], "vi": r["VIvega"], "x": r["x"], "y": r["y"]})


# 1a. his own output
r2 = rawframe(M2)
z2 = martin_zero_point(r2, M2FID, 18.85)
f = pd.read_csv(f"{MARTIN}/M15_fitstars_ZPCR.dat", index_col=0)
d, j = cKDTree(r2[["x", "y"]].values).query(f[["x", "y"]].values)
diff = f.vi.values - (r2.vi.values[j] - z2[j])
print(f"1a. Martin's M15_fitstars_ZPCR.dat ({len(f)} stars, all matched within {d.max():.3g} px): "
      f"his vi - ours: max |.| {np.abs(diff).max():.2e}, median {np.median(diff):.2e}")

# 1b. the notebook's loop on our M15 data
r = rawframe(f"{MARTIN}/M15.dat")
s_m = json.load(open(f"{B}/m15_martinzp/setup.json"))
z = martin_zero_point(r, s_m["fiducial"], s_m["sgb_v"], s_m["zp_k"], s_m["zp_martin_dmag"])
use = np.flatnonzero((r.v < s_m["sgb_v"] + 3.5) & (r.v > s_m["sgb_v"] - 3.5))
fid = pd.read_csv(s_m["fiducial"])
ff = CubicSpline(fid["v"], fid["vi"])
xy = r[["x", "y"]].values[use]
dc = r.vi.values[use] - ff(r.v.values[use])
sub = np.random.default_rng(0).choice(len(use), 300, replace=False)
bf = np.array([np.median(dc[np.argpartition(np.sum((xy[q] - xy) ** 2, axis=1), 100)[:100]]) for q in sub])
print(f"1b. notebook loop vs martin_zero_point, 300 stars of our M15.dat: max |diff| {np.abs(bf - z[use][sub]).max():.2e}")

# 2. production reproduced
lim = {int(k): v for k, v in json.load(open("/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json")).items()}
g = pd.read_csv(f"{B}/m15_final_refined/combined_2d_ks_results.csv")
g = g[g.age_gyr < g.simulation_id.map(lim) - 1e-6]
row = g.loc[g.ks_statistic.idxmin()]
run, age, mu, e = int(row.simulation_id), float(row.age_gyr), float(row.distance_modulus), float(row.extinction_ebv)


def ks_at(setup_path, override=None):
    d = {k: v for k, v in json.load(open(setup_path)).items() if k in Setup.__dataclass_fields__}
    d.update(override or {})
    S = Setup(**d)
    t0 = time.time()
    ctx = context(S)
    _, make = ctx.model(run, age)
    cat = make(np.random.RandomState(S.crn_seed), mu, e)
    out = ctx.stage3(cat)
    return ctx, out, time.time() - t0


ov = json.load(open(f"{B}/m15_final/refine/setup_override.json"))
ctx_p, out_p, t = ks_at(f"{B}/m15_final_refined/setup.json", ov)
print(f"2. production cell run {run}, {age} Gyr, mu {mu}, E {e}: file KS {row.ks_statistic:.6f}, "
      f"this code (zp_method ours) {out_p['ks_statistic']:.6f}, n_obs {out_p['n_observed']} ({t:.0f} s)")

# 3. Martin's map
ctx_m, out_m, t = ks_at(f"{B}/m15_martinzp/setup.json", ov)
zp_p, zp_m = ctx_p.zp_map, ctx_m.zp_map
win = (r.v > s_m["sgb_v"] - 2) & (r.v < s_m["sgb_v"] + 2)
for name, m in (("window 16.64-20.64, all stars", win.values),):
    a, b = zp_p[m], zp_m[m]
    print(f"3. colour zero point over {name} ({m.sum()} stars): ours rms {a.std():.4f} (median {np.median(a):+.4f}), "
          f"Martin rms {b.std():.4f} (median {np.median(b):+.4f}); correlation {np.corrcoef(a, b)[0, 1]:.3f}; "
          f"Martin - ours rms {np.std(b - a):.4f}, median {np.median(b - a):+.4f}")
print(f"   F606W: ours moves it by 0.57 x the colour correction, Martin keeps it fixed")
print(f"   fit stars in the window: production {ctx_p.n_target}, Martin zp {ctx_m.n_target}")
print(f"   same cell with Martin's zero point: KS {out_m['ks_statistic']:.6f} ({t:.0f} s)")
np.savez(f"{B}/m15_martinzp/analysis/zp_maps.npz", x=r.x.values, y=r.y.values, v=r.v.values, vi=r.vi.values,
         ours=zp_p, martin=zp_m)
