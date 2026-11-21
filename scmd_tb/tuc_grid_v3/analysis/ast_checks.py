#!/usr/bin/env python3
"""Diagnostics of the 47 Tuc artificial-star error model (v3 settings).

1. Shuffle test: inside each synthetic star's own k = 30 neighbourhood, take the
   F606W error from one recovered neighbour and the F814W error from another.
   That keeps the conditioning and breaks only the pairing, so it separates the
   effect of pairing from the effect of MakeTestCMD's coarse 10 x 12 bins.
2. AST-only colour map: output - input colour of recovered artificial stars on
   the upper MS, median of the k nearest on the detector. Should be flat if the
   real-star colour map (colour_zero_point) is a PSF/calibration effect.
3. CTE: trend of (input - output) with detector y for bright, uncrowded ASTs.
4. Observed catalogue: are there signs of upstream quality cuts the ASTs lack?
"""
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_v3/code")
from isochronetoolbox.scmd.ast import ArtificialStars, centre  # noqa: E402

M = "/dartfs-hpc/scratch/f007058/Files_from_Martin/47Tuc_data"
ZPV, ZPI = 30.5742, 29.6777
rob = lambda q: 1.4826 * np.median(np.abs(q - np.median(q)))
rs = np.random.RandomState(20261031)

a = pd.read_csv(f"{M}/47Tucartstars.dat")
fit = pd.read_csv(f"{M}/47Tuc_fitstars.dat")
win = {"F606W": (15.0, 19.1), "F814W": (14.0, 18.6)}
A = ArtificialStars(a, fit, {"F606W": ZPV, "F814W": ZPI}, win, {"F606W": 16.48, "F814W": 15.55},
                    k=30, smooth=True, metric="crowding", radial="observed")
ne = A.nearest
x, y = a.InputX.values, a.InputY.values
xc, yc = centre(x, y)
r = np.hypot(x - xc, y - yc)
V = a.InputF606W.values + ZPV
q = np.flatnonzero((V > 15.0) & (V < 19.0))

# ---------------------------------------------------------------- 1. shuffle test
nb = ne.nbr[q]                                    # each query's 30 neighbours (query = an AST input)
ok = ne.rec[nb]
has = ok.sum(axis=1) >= 2
nb, ok, qq = nb[has], ok[has], q[has]


def pick(mask):
    s = np.where(mask, rs.random_sample(mask.shape), -1.0)
    return np.take_along_axis(nb, s.argmax(axis=1)[:, None], axis=1)[:, 0]


p1, p2 = pick(ok), pick(ok)
col = {"paired (nearest)": ne.ev[p1] - ne.ei[p1],
       "unpaired, same neighbourhood": ne.ev[p1] - ne.ei[p2]}
okm, evm, eim = ne.draw(rs, r[qq], V[qq] - ZPV, a.InputF814W.values[qq])
col["smooth (v3 default)"] = (evm - eim)[okm]
rad_smooth = r[qq][okm]
okv, ev = A.band["F606W"].draw(rs, r[qq], V[qq] - ZPV)
oki, ei = A.band["F814W"].draw(rs, r[qq], a.InputF814W.values[qq])
col["MakeTestCMD bins, independent"] = (ev - ei)[okv & oki]
rad_bins = r[qq][okv & oki]
radius = {"paired (nearest)": r[qq], "unpaired, same neighbourhood": r[qq],
          "smooth (v3 default)": rad_smooth, "MakeTestCMD bins, independent": rad_bins}
print("1. colour error (dV - dI) of recovered stars, F606W 15-19: robust sigma [mmag] / tail fraction |dcol| > 0.05")
rings = [(0, 300), (300, 600), (600, 1000), (1000, 1600), (1600, 4000)]
print(f"   {'':32s}" + "".join(f"{f'{lo}-{hi} px':>17s}" for lo, hi in rings) + f"{'all':>17s}")
for lab, c in col.items():
    cells = []
    for lo, hi in rings + [(0, 1e9)]:
        s = (radius[lab] >= lo) & (radius[lab] < hi)
        cc = c[s]
        cells.append(f"{1000 * rob(cc):7.1f} / {np.mean(np.abs(cc - np.median(cc)) > 0.05):.3f}")
    print(f"   {lab:32s}" + "".join(f"{v:>17s}" for v in cells))

# ---------------------------------------------------------------- 2. AST-only colour map
dv = (a.InputF606W - a.OutputF606W).values
di = (a.InputF814W - a.OutputF814W).values
rec = ne.rec
sel = rec & (V > 17.6) & (V < 18.9)
dcol = -(dv - di)[sel]                            # output - input colour
pts = np.column_stack([x[sel], y[sel]])
_, idx = cKDTree(pts).query(pts, k=101)
amap = np.median(dcol[idx[:, 1:]], axis=1)        # k = 100 nearest, self excluded
null = np.median(rs.permutation(dcol)[idx[:, 1:]], axis=1)   # same, positions scrambled
print("\n2. AST-only colour map, recovered ASTs F606W 17.6-18.9, median of 100 nearest on the detector [mmag]")
for lab, mm in (("map", amap), ("positions scrambled (noise)", null)):
    p = np.percentile(mm, [1, 50, 99]) * 1000
    print(f"   {lab:28s} std {1000 * mm.std():5.2f}   1-99% {p[0]:+6.1f} to {p[2]:+6.1f}   median {p[1]:+5.1f}")
outer = np.hypot(pts[:, 0] - xc, pts[:, 1] - yc) > 1000
print(f"   map std outside 1000 px: {1000 * amap[outer].std():.2f}  (noise {1000 * null[outer].std():.2f})")

# ---------------------------------------------------------------- 3. CTE
print("\n3. CTE: slope of (input - output) with detector y, recovered ASTs F606W < 17.5, r > 1200 px [mmag per 1000 px]")
b = rec & (V < 17.5) & (r > 1200)
for lab, d in (("F606W", dv), ("F814W", di), ("colour", dv - di)):
    for coord, cname in ((y, "y"), (x, "x")):
        A_ = np.column_stack([np.ones(b.sum()), coord[b] / 1000.0])
        coef, *_ = np.linalg.lstsq(A_, d[b], rcond=None)
        res = d[b] - A_ @ coef
        se = rob(res) / np.sqrt(b.sum()) / (coord[b] / 1000.0).std()
        print(f"   {lab:7s} vs {cname}: {1000 * coef[1]:+6.2f} +- {1000 * se:4.2f}", end="")
    print()
print(f"   ({b.sum()} stars; median (input - output): F606W {1000 * np.median(dv[b]):+.1f}, F814W {1000 * np.median(di[b]):+.1f} mmag)")

# ---------------------------------------------------------------- 4. catalogue quality columns
obs = pd.read_csv(f"{M}/47Tuc.dat", sep=r"\s+")
w = obs[(obs.Vvega > 15.0) & (obs.Vvega < 19.0)]
print(f"\n4. observed catalogue 47Tuc.dat, F606W 15-19: {len(w)} stars (of {len(obs)})")
for c in ("qfitV", "qfitI", "othv", "othi", "Nv", "Ni", "xsig", "ysig", "err", "err.2"):
    if c in w:
        p = np.percentile(w[c], [0, 1, 50, 99, 100])
        print(f"   {c:6s} min {p[0]:9.4f}  1% {p[1]:9.4f}  median {p[2]:9.4f}  99% {p[3]:9.4f}  max {p[4]:9.4f}")
