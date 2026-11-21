#!/usr/bin/env python3
"""Error terms the artificial stars cannot contain, measured from Treasury - HUGS.

HUGS re-reduced the same GO-10775 F606W/F814W images, so in the star-by-star
difference astrophysical structure (reddening, populations) cancels, and shared photon
noise and blends largely cancel; what is left is reduction-specific error.

  map       median (dF606W, dF814W) of the k = 50 nearest matched stars on the detector
            (self excluded): the spatial part. Compared with the CMD-residual colour map
            of v2 (colour offset from the ridge minus the AST-predicted offset), it shows
            how much of that map is reduction-specific.
  per-star  (dF606W, dF814W) about the map, by F606W and radius: the random part, kept
            as pairs. It also contains HUGS's own errors, so it is an upper estimate
            (/sqrt 2 if the two reductions' own errors are equal and independent).

Writes hugs_terms.npz (matched stars: x, y, F606W, radius, map and residual pairs).
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

F = os.path.dirname(os.path.abspath(__file__))
B = "/dartfs-hpc/scratch/f007058/scmd_tb"
sys.path.insert(0, os.path.join(F, "code"))
from isochronetoolbox.scmd.ast import centre  # noqa: E402

rob = lambda q: 1.4826 * np.median(np.abs(q - np.median(q)))
T = pd.read_csv("/dartfs-hpc/scratch/f007058/Files_from_Martin/47Tuc_data/47Tuc.dat", sep=r"\s+")
H = pd.read_csv(f"{B}/hugs/47Tuc_hugs_meth2.dat", sep=r"\s+")
zpmap = np.load(f"{B}/tuc_grid_v2/analysis/zp_map.npy")          # one value per Treasury star
assert len(zpmap) == len(T)

d, j = cKDTree(H[["x", "y"]].values).query(T[["x", "y"]].values)
V, I, C = T.Vvega.values, T.Ivega.values, T.VIvega.values
dv = V - H.Vvega.values[j]
di = I - H.Ivega.values[j]
m = (d < 0.5) & (np.abs(dv) < 0.1) & (np.abs(di) < 0.1) & (V > 15.0) & (V < 19.1)
# stars on the sequence: |colour - running median| < 0.08 (drops blends, binaries, field)
edges = np.arange(15.0, 19.2, 0.1)
kb = np.clip(np.digitize(V, edges) - 1, 0, len(edges) - 2)
med = np.array([np.median(C[m & (kb == b)]) for b in range(len(edges) - 1)])
m &= np.abs(C - med[kb]) < 0.08
idx = np.flatnonzero(m)
xc, yc = centre(T.x.values, T.y.values)
r = np.hypot(T.x.values - xc, T.y.values - yc)
a, b = dv[idx] - np.median(dv[idx]), di[idx] - np.median(di[idx])
xy = np.column_stack([T.x.values[idx], T.y.values[idx]])
print(f"matched sequence stars, F606W 15-19.1: {len(idx)} of {int(((V > 15) & (V < 19.1)).sum())} Treasury stars")

K = 50
_, nb = cKDTree(xy).query(xy, k=K + 1)
nb = nb[:, 1:]
mv, mi = np.median(a[nb], axis=1), np.median(b[nb], axis=1)
rs = np.random.RandomState(1)
pa, pb = rs.permutation(a), rs.permutation(b)
nv, ni = np.median(pa[nb], axis=1), np.median(pb[nb], axis=1)
print("\n1. Treasury - HUGS map (median of 50 nearest) [mmag]:  std  /  1-99%  /  noise std (positions scrambled)")
for lab, q, n in (("F606W", mv, nv), ("F814W", mi, ni), ("colour", mv - mi, nv - ni)):
    p = np.percentile(q, [1, 99]) * 1000
    print(f"   {lab:7s} {1000 * q.std():5.2f}   {p[0]:+6.1f}..{p[1]:+6.1f}   noise {1000 * n.std():4.2f}")
zc = zpmap[idx]
cm = mv - mi
rho = np.corrcoef(cm, zc)[0, 1]
slope = np.polyfit(cm, zc, 1)[0]
print(f"\n2. CMD-residual colour map (v2) at the same stars: std {1000 * zc.std():.2f} mmag")
print(f"   correlation with the Treasury - HUGS colour map: {rho:+.2f};  slope (CMD map per T-H map) {slope:+.2f}")
print(f"   CMD map variance explained by the T-H colour map: {rho ** 2:.2f}  ->  reduction-specific part ~{1000 * abs(rho) * zc.std():.1f} mmag,"
      f" the rest ~{1000 * np.sqrt(max(zc.var() * (1 - rho ** 2), 0)):.1f} mmag")
for lo_r, hi_r in ((0, 600), (600, 1200), (1200, 1e9)):
    s = (r[idx] >= lo_r) & (r[idx] < hi_r)
    print(f"   r {lo_r:4.0f}-{min(hi_r, 9999):4.0f}: T-H colour map std {1000 * cm[s].std():4.2f}, CMD map std {1000 * zc[s].std():4.2f},"
          f" correlation {np.corrcoef(cm[s], zc[s])[0, 1]:+.2f}")

ra, rb = a - mv, b - mi
print("\n3. per-star (dF606W, dF814W) about the map, robust sigma [mmag] and correlation r (|z| < 5)")
print("   zone        F606W      n   sigV  sigI  sig_col     r")
rows = []
for lo_r, hi_r in ((0, 600), (600, 1200), (1200, 1e9)):
    for lo in np.arange(15.0, 19.01, 0.5):
        s = (r[idx] >= lo_r) & (r[idx] < hi_r) & (V[idx] >= lo) & (V[idx] < lo + 0.5)
        if s.sum() < 50:
            continue
        x1, x2 = ra[s], rb[s]
        core = (np.abs(x1 - np.median(x1)) < 5 * rob(x1)) & (np.abs(x2 - np.median(x2)) < 5 * rob(x2))
        rr = np.corrcoef(x1[core], x2[core])[0, 1]
        rows.append((lo_r, lo, s.sum(), rob(x1), rob(x2), rob(x1 - x2), rr))
        print(f"   r {lo_r:4.0f}-{min(hi_r, 9999):4.0f}  {lo:4.1f}-{lo + 0.5:4.1f}  {s.sum():6d}  {1000 * rob(x1):4.1f}  {1000 * rob(x2):4.1f}  {1000 * rob(x1 - x2):6.1f}  {rr:+.2f}")
np.savez(os.path.join(F, "hugs_terms.npz"), x=xy[:, 0], y=xy[:, 1], v=V[idx], r=r[idx], map_v=mv, map_i=mi,
         res_v=ra, res_i=rb, xc=xc, yc=yc)
print(f"\n{os.path.join(F, 'hugs_terms.npz')}")
