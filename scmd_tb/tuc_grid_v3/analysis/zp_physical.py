#!/usr/bin/env python3
"""Is the Treasury-catalogue colour pattern physical?

1. HUGS (independent reduction of the same F606W/F814W images): cross-match star by
   star, map the colour difference Treasury - HUGS on the upper main sequence, and
   compare it with the zero-point map the v2 grid removed.  Also each catalogue's
   own colour map.  A noise floor from shuffled positions.
2. Reddening vector: differential reddening moves a star in F606W by ~2.8x its
   colour shift (after the subgiant branch's tilt is taken out); an F814W zero-point error moves it not at all, an F606W
   one by 1x.  Regress the SGB magnitude map on the colour map.
"""
import json
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
B = "/dartfs-hpc/scratch/f007058/scmd_tb"
rng = np.random.default_rng(1)
acs = pd.read_csv(json.load(open(f"{B}/tuc_grid_v2/setup.json"))["obs_raw"], sep=r"\s+")
zp = np.load(f"{B}/tuc_grid_v2/analysis/zp_map.npy")
hug = pd.read_csv(f"{B}/hugs/47Tuc_hugs_meth2.dat", sep=r"\s+")
ax_, ay_, av, avi = acs.x.values, acs.y.values, acs.Vvega.values, acs.VIvega.values
hx, hy, hv, hvi = hug.x.values, hug.y.values, hug.Vvega.values, hug.VIvega.values
K = 50


def local_median(x, y, val, k=K, at=None):
    """median of val over the k nearest (x, y) points (self excluded when at is None)"""
    t = cKDTree(np.column_stack([x, y]))
    if at is None:
        _, i = t.query(np.column_stack([x, y]), k=k + 1); i = i[:, 1:]
    else:
        _, i = t.query(np.column_stack(at), k=k)
    return np.median(val[i], axis=1)


def ms_resid(v, vi, sel):
    """colour minus a clipped polynomial ridge colour(V) on the upper MS"""
    m = sel.copy()
    for _ in range(5):
        c = np.polyfit(v[m], vi[m], 4)
        r = vi - np.polyval(c, v)
        s = 1.4826 * np.median(np.abs(r[m] - np.median(r[m])))
        m = sel & (np.abs(r) < 3 * s)
    return r, m


# ---- 1. cross-match ----------------------------------------------------------------
t = cKDTree(np.column_stack([hx, hy]))
d, j = t.query(np.column_stack([ax_, ay_]))
ok = (d < 0.5) & (np.abs(av - hv[j]) < 0.1)
print(f"matched {ok.sum()} of {len(av)} Treasury stars (0.5 px, 0.1 mag); median sep {np.median(d[ok]):.3f} px")
ums = ok & (av > 17.6) & (av < 18.9)
dcol = avi - hvi[j]; dv = av - hv[j]
dcol_c = dcol - np.median(dcol[ums]); dv_c = dv - np.median(dv[ums])
print(f"upper-MS matched {ums.sum()}: colour diff median {np.median(dcol[ums]):+.4f}, "
      f"per-star scatter {1.4826*np.median(np.abs(dcol_c[ums])):.4f}")
idx = np.flatnonzero(ums)
x, y = ax_[idx], ay_[idx]
diff_map = local_median(x, y, dcol_c[idx])
diff_noise = local_median(x, y, rng.permutation(dcol_c[idx]))
ra, ma = ms_resid(av, avi, (av > 17.6) & (av < 18.9))
rh, mh = ms_resid(hv, hvi, (hv > 17.6) & (hv < 18.9))
both = ums & ma & mh[j]
ib = np.flatnonzero(both)
acs_map = local_median(ax_[ib], ay_[ib], ra[ib])
hug_map = local_median(ax_[ib], ay_[ib], rh[j][ib])
acs_noise = local_median(ax_[ib], ay_[ib], rng.permutation(ra[ib]))
zp_at = zp[idx]


def corr(a, b):
    return np.corrcoef(a, b)[0, 1]


def slope(xv, yv):
    return np.polyfit(xv, yv, 1)[0]


print("\n1. maps over the matched upper MS (std in mmag; noise = same median of shuffled positions)")
print(f"   Treasury own colour map      std {acs_map.std()*1e3:5.1f}   (noise {acs_noise.std()*1e3:.1f})")
print(f"   HUGS own colour map          std {hug_map.std()*1e3:5.1f}")
print(f"   Treasury - HUGS colour map   std {diff_map.std()*1e3:5.1f}   (noise {diff_noise.std()*1e3:.1f})")
print(f"   v2 zero-point map (these stars) std {zp_at.std()*1e3:5.1f}")
print(f"   corr(Treasury - HUGS, zp map) = {corr(diff_map, zp_at):.2f}, slope {slope(zp_at, diff_map):.2f}")
print(f"   corr(HUGS own map, Treasury own map) = {corr(hug_map, acs_map):.2f}")
print(f"   corr(HUGS own map, zp map) = {corr(hug_map, zp[ib]):.2f}")

# ---- 2. reddening vector on the subgiant branch ----------------------------------------
sgb = (avi > 0.57) & (avi < 0.66) & (av > 16.75) & (av < 17.25)
m = sgb.copy()
for _ in range(5):
    c = np.polyfit(avi[m], av[m], 2); r = av - np.polyval(c, avi)
    s = 1.4826 * np.median(np.abs(r[m] - np.median(r[m]))); m = sgb & (np.abs(r) < 3 * s)
isg = np.flatnonzero(m)
ks = 60
# a colour shift dc and an F606W shift dV move a star's SGB residual by dV - s*dc,
# s = dV_ridge/dcolour at its colour (the SGB is tilted); with dV = beta*dc the
# tilt-corrected residual r + s*dc = beta*dc
srid = np.polyval(np.polyder(c), avi[isg])
tq = cKDTree(np.column_stack([ax_[isg], ay_[isg]]))
_, nb = tq.query(np.column_stack([ax_[isg], ay_[isg]]), k=ks + 1); nb = nb[:, 1:]
col_at_sgb = zp[isg]
sgb_map = np.median(r[isg][nb], axis=1) + np.median(srid[nb], axis=1) * col_at_sgb
beta = np.sum(sgb_map * col_at_sgb) / np.sum(col_at_sgb ** 2)
print(f"\n2. subgiant branch: {len(isg)} stars, ridge tilt dV/dcolour median {np.median(srid):.2f}")
print(f"   beta = dF606W / dcolour = {beta:.2f}, corr {corr(col_at_sgb, sgb_map):.2f}")
print("   expected: reddening ~ +2.8 (A606/E(606-814)), F606W zero point +1, F814W zero point 0")

fig, axs = plt.subplots(1, 4, figsize=(20, 5.2))
kw = dict(s=1.2, cmap="RdBu_r", vmin=-12, vmax=12, rasterized=True)
for a, (val, xx, yy, tt) in zip(axs[:3], [(zp_at, x, y, "Treasury: zero-point map removed in v2"),
                                          (hug_map, ax_[ib], ay_[ib], "HUGS (same images): its own colour map"),
                                          (diff_map, x, y, "Treasury - HUGS, star by star")]):
    sc = a.scatter(xx, yy, c=val * 1e3, **kw); a.set_aspect("equal"); a.set_title(tt, fontsize=10)
    a.set_xlabel("x (px)")
axs[0].set_ylabel("y (px)")
fig.colorbar(sc, ax=axs[:3].tolist(), shrink=0.8, pad=0.01, label="colour offset (mmag)")
a = axs[3]
a.scatter(col_at_sgb * 1e3, sgb_map * 1e3, s=2, c="k", alpha=0.4, rasterized=True)
xx = np.linspace(-15, 15, 2)
for sl, lab, cc in [(2.8, "differential reddening (2.8)", "#c0392b"), (1.0, "F606W zero point (1)", "#e9a23b"), (0.0, "F814W zero point (0)", "#2f6db5")]:
    a.plot(xx, sl * xx, color=cc, lw=1.6, label=lab)
a.set_xlabel("colour offset at the star (mmag)"); a.set_ylabel("subgiant-branch F606W offset, tilt-corrected (mmag)")
a.set_title(f"reddening-vector test: fitted slope {beta:.2f}", fontsize=10)
a.set_xlim(-18, 18); a.set_ylim(-40, 40); a.legend(fontsize=8)
fig.savefig(f"{B}/tuc_grid_v3/analysis/zp_physical.png", dpi=110, bbox_inches="tight")
print(f"{B}/tuc_grid_v3/analysis/zp_physical.png")
