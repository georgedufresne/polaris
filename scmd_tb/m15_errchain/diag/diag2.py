#!/usr/bin/env python3
"""Is M15's spatial colour pattern differential reddening or a zero-point pattern?
Per spatial cell, fit a shift (d colour, d F606W) of the observed ridge to the cell's stars: the steep
upper MS fixes mostly the colour shift, the flat SGB mostly the magnitude shift.  Across cells the shifts
line up along dF606W/dcolour = A_F606W / E(F606W-F814W) ~ 2.8 for reddening, 0 for an F814W zero point,
1 for an F606W zero point."""
import json
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

E = "/dartfs-hpc/scratch/f007058/scmd_tb/m15_errchain"
sys.path.insert(0, f"{E}/code")
from isochronetoolbox.scmd import Setup, context             # noqa: E402

S = Setup(**{k: v for k, v in json.load(open(f"{E}/setup.json")).items() if k in Setup.__dataclass_fields__})
C = context(S)
fit = C.ast._args[0]
xc, yc = C.ast.centre_xy
v, c, x, y = fit.v.values, (fit.v - fit.i).values, fit.x.values, fit.y.values
r = np.hypot(x - xc, y - yc)
rid = C.rcut.ridge
# upper MS: colour as a function of F606W; SGB: F606W as a function of colour
m1 = rid.mag > 19.3
o1 = np.argsort(rid.mag[m1])
ms_c = lambda V: np.interp(V, rid.mag[m1][o1], rid.color[m1][o1])
m2 = (rid.mag > 18.2) & (rid.mag < 19.0) & (rid.color > 0.52) & (rid.color < 0.66)
o2 = np.argsort(rid.color[m2])
sgb_v = lambda col: np.interp(col, rid.color[m2][o2], rid.mag[m2][o2])
MS = (v > 19.4) & (v < 20.6)
SG = (c > 0.54) & (c < 0.63) & (v > 18.2) & (v < 19.0) & (np.abs(v - sgb_v(c)) < 0.15)
h = 0.01
s_ms = np.median((ms_c(v[MS] + h) - ms_c(v[MS] - h)) / (2 * h))     # d colour / d F606W on the MS
s_sg = np.median((sgb_v(c[SG] + h) - sgb_v(c[SG] - h)) / (2 * h))   # d F606W / d colour on the SGB
print(f"MS stars {MS.sum()}, SGB stars {SG.sum()}; ridge slopes: MS dcol/dV {s_ms:.4f}, SGB dV/dcol {s_sg:.2f}")


def rsig(a):
    return 1.4826 * np.median(np.abs(a - np.median(a)))


rows = []
NB = 7
far = r > 400                                   # leave the crowded core out
xe = np.quantile(x[far], np.linspace(0, 1, NB + 1))
ye = np.quantile(y[far], np.linspace(0, 1, NB + 1))
for p in range(NB):
    for q in range(NB):
        cell = far & (x >= xe[p]) & (x <= xe[p + 1]) & (y >= ye[q]) & (y <= ye[q + 1])
        a, b = cell & MS, cell & SG
        if a.sum() < 200 or b.sum() < 25:
            continue
        ra = c[a] - ms_c(v[a])                   # = dc - s_ms dV
        rb = v[b] - sgb_v(c[b])                  # = dV - s_sg dc
        m_a, e_a = np.median(ra), 1.25 * rsig(ra) / np.sqrt(a.sum())
        m_b, e_b = np.median(rb), 1.25 * rsig(rb) / np.sqrt(b.sum())
        A = np.array([[1.0, -s_ms], [-s_sg, 1.0]])
        Ai = np.linalg.inv(A)
        dc, dv = Ai @ np.array([m_a, m_b])
        cov = Ai @ np.diag([e_a ** 2, e_b ** 2]) @ Ai.T
        rows.append((np.mean(x[cell]), np.mean(y[cell]), dc, dv, np.sqrt(cov[0, 0]), np.sqrt(cov[1, 1]),
                     cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1]), a.sum(), b.sum()))
t = np.array(rows)
dc, dv, ec, ev = t[:, 2] - np.median(t[:, 2]), t[:, 3] - np.median(t[:, 3]), t[:, 4], t[:, 5]
print(f"{len(t)} cells (r > 400 px): colour shifts rms {np.std(dc):.4f} (noise {np.median(ec):.4f}), "
      f"F606W shifts rms {np.std(dv):.4f} (noise {np.median(ev):.4f})")


def chi2(k):
    """chi2 of dv = k dc + const, errors in both (with their correlation)."""
    rho = t[:, 6]
    var = ev ** 2 + k ** 2 * ec ** 2 - 2 * k * rho * ec * ev
    resid = dv - k * dc
    w = 1 / var
    z = np.sum(w * resid) / np.sum(w)
    return np.sum(w * (resid - z) ** 2)


ks = np.linspace(-3, 8, 1101)
cs = np.array([chi2(k) for k in ks])
kb = ks[cs.argmin()]
inside = ks[cs <= cs.min() + 1]
print(f"best slope dF606W/dcolour = {kb:.2f} (1 sigma {inside.min():.2f} to {inside.max():.2f}); "
      f"chi2/dof {cs.min() / (len(t) - 2):.2f}")
for name, k in (("reddening (2.8)", 2.8), ("F606W zero point (1)", 1.0), ("F814W zero point (0)", 0.0)):
    print(f"  {name:22s}: chi2 {chi2(k):7.1f}, delta chi2 vs best {chi2(k) - cs.min():7.1f}")

fig, axs = plt.subplots(1, 2, figsize=(13, 5.5), constrained_layout=True)
ax = axs[0]
ax.errorbar(dc, dv, xerr=ec, yerr=ev, fmt="o", ms=4, color="#2166ac", ecolor="0.7", lw=0.8)
xx = np.array([dc.min() - 0.003, dc.max() + 0.003])
for k, ls, lab in ((kb, "-", f"best fit {kb:.2f}"), (2.8, "--", "reddening 2.8"), (1.0, ":", "F606W zero point 1"),
                   (0.0, "-.", "F814W zero point 0")):
    ax.plot(xx, k * xx, ls, color="0.2" if ls == "-" else "0.5", lw=1.2, label=lab)
ax.set_xlabel("cell colour shift (mag)"); ax.set_ylabel("cell F606W shift (mag)"); ax.legend(fontsize=8)
ax.set_title("per-cell shift of the CMD (MS + SGB), r > 400 px", fontsize=10)
ax = axs[1]
lim = np.max(np.abs(dc))
sc = ax.scatter(t[:, 0], t[:, 1], c=dc, cmap="RdBu_r", vmin=-lim, vmax=lim, s=260, marker="s")
for xi, yi, dvi in zip(t[:, 0], t[:, 1], dv):
    ax.text(xi, yi, f"{1000 * dvi:+.0f}", ha="center", va="center", fontsize=7)
ax.plot(xc, yc, "k+", ms=12)
fig.colorbar(sc, ax=ax, label="colour shift (mag); numbers: F606W shift (mmag)")
ax.set_aspect("equal"); ax.set_xlabel("x (px)"); ax.set_ylabel("y (px)")
fig.suptitle("M15: differential reddening or zero-point pattern?")
fig.savefig(f"{E}/diag/reddening_vs_zeropoint.png", dpi=110)
print(f"{E}/diag/reddening_vs_zeropoint.png")
