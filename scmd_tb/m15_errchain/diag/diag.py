#!/usr/bin/env python3
"""Error-model diagnostics for the M15 sCMD: radial sample, completeness floor, neighbourhood reach,
sign, zero points, radial closure, and where the observed CMD is wider than the model."""
import json
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

E = "/dartfs-hpc/scratch/f007058/scmd_tb/m15_errchain"
sys.path.insert(0, f"{E}/code")
from isochronetoolbox.scmd import Setup, context             # noqa: E402
from isochronetoolbox.scmd.ast import centre                 # noqa: E402

OUT = f"{E}/diag"
BASE = json.load(open(f"{E}/setup.json"))
RUN, AGE, MU, EV = 10286, 14.1, 15.48, 0.08                  # the grid's best cell
RB = np.array([0, 100, 200, 400, 700, 1200, 2000, 4500])     # radius bins (px)


def ctx(**kw):
    d = {k: v for k, v in BASE.items() if k in Setup.__dataclass_fields__}
    d.update(kw)
    return context(Setup(**d))


def rsig(x):
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) > 5 else np.nan


def rbin(r):
    return np.clip(np.digitize(r, RB) - 1, 0, len(RB) - 2)


def rlab(k):
    return f"{RB[k]:4d}-{RB[k + 1]:4d}"


C = ctx()                                         # production: radii from the fit stars (after the cut)
R = ctx(radial_source="raw")                      # radii from the loose-cut photometry
S = C.s
zp = {"v": S.zp_v, "i": S.zp_i}
lo, hi = C.v_window
fit = C.ast._args[0]
xc, yc = C.ast.centre_xy
fit_r = np.hypot(fit.x.values - xc, fit.y.values - yc)
loose = R.radial_stars.copy()
loose["r"] = np.hypot(loose.x.values - xc, loose.y.values - yc)
print(f"window F606W {lo:.2f}-{hi:.2f}; fit stars {len(fit)}; loose-cut stars {len(loose)} "
      f"({((loose.v > lo) & (loose.v < hi)).sum()} in the window); centre ({xc:.1f}, {yc:.1f})")

# ---- 1. the ridge cut's efficiency s(r, m): fraction of loose-cut stars that survive it -----------
inwin = loose[(loose.v > lo) & (loose.v < hi)].reset_index(drop=True)
inwin = inwin.merge(fit[["x", "y"]].drop_duplicates().assign(kept=True), on=["x", "y"], how="left")
inwin["kept"] = inwin["kept"].fillna(False).astype(bool)
print(f"fit stars matched in the loose-cut sample: {inwin.kept.sum()} of {((fit.v > lo) & (fit.v < hi)).sum()}")
mb = [(lo, 17.6, "upper RGB"), (17.6, 18.5, "lower RGB"), (18.5, 19.0, "SGB"), (19.0, 19.7, "TO"), (19.7, hi, "MS")]
print("\n[1] ridge-cut efficiency s(r, m) = kept / loose-cut stars")
print("  r (px)     " + "".join(f"{n:>11s}" for _, _, n in mb))
for k in range(len(RB) - 1):
    row = []
    for a, b, _ in mb:
        sel = (rbin(inwin.r.values) == k) & (inwin.v > a) & (inwin.v <= b)
        row.append(f"{inwin.kept[sel].mean():6.3f}({sel.sum():5d})" if sel.sum() > 20 else f"{'-':>11s}")
    print(f"  {rlab(k)}  " + "".join(row))
print(f"  all: {inwin.kept.mean():.3f}")

# ---- 2. completeness floor ----------------------------------------------------------------------------
N = C.ast.nearest
for name, df, r in (("fit stars", fit, fit_r), ("loose-cut stars", loose, loose.r.values)):
    c = N.completeness(df.v.values - zp["v"], df.i.values - zp["i"], r)
    print(f"\n[2] {name}: completeness (100 nearest ASTs) < 0.1 for {np.mean(c < 0.1):.4f}, < 0.2 for {np.mean(c < 0.2):.4f}")
    print("  r (px)     median c   frac c<0.1   n")
    for k in range(len(RB) - 1):
        sel = rbin(r) == k
        if sel.sum():
            print(f"  {rlab(k)}   {np.median(c[sel]):7.3f}   {np.mean(c[sel] < 0.1):9.4f}  {sel.sum():6d}")

# ---- 3. neighbourhood reach; the nearest-AST shortcut vs an exact query --------------------------------
a = C.ast._full
ax_, ay_ = centre(a["InputX"].values, a["InputY"].values)
ar = np.hypot(a["InputX"].values - ax_, a["InputY"].values - ay_)
avin = a["InputF606W"].values.astype(float)
rs = np.random.RandomState(1)
smp = inwin.iloc[rs.choice(len(inwin), min(20000, len(inwin)), replace=False)]
co = N._coords(smp.v.values - zp["v"], smp.i.values - zp["i"], smp.r.values)
_, near = N.tree.query(co)
nb_s = N.nbr[near]
_, nb_e = N.tree.query(co, k=N.k)
vst = smp.v.values - zp["v"]
print(f"\n[3] neighbourhoods of {len(smp)} in-window stars (k = {N.k})")
print("  r (px)    ASTs/obs  |dV| med  dV span  r_nb/r p10-p90   overlap  sig(dcol) shortcut / exact  n")
rows3 = []
for k in range(len(RB) - 1):
    sel = rbin(smp.r.values) == k
    if sel.sum() < 30:
        continue
    dv = np.abs(avin[nb_e[sel]] - vst[sel, None])
    rr = ar[nb_e[sel]] / np.maximum(smp.r.values[sel, None], 1.0)
    ov = np.mean([len(np.intersect1d(p, q)) / N.k for p, q in zip(nb_s[sel], nb_e[sel])])
    sig = []
    for nb in (nb_s[sel], nb_e[sel]):
        ok = N.rec[nb]
        dcol = np.where(ok, N.ev[nb] - N.ei[nb], np.nan)
        sig.append(np.nanmedian([rsig(d) for d in dcol]))
    dens = (rbin(ar) == k).sum() / max((rbin(inwin.r.values) == k).sum(), 1)
    rows3.append((RB[k], np.median(dv), np.median(np.ptp(avin[nb_e[sel]], axis=1)), np.percentile(rr, 10),
                  np.percentile(rr, 90), ov, sig[0], sig[1]))
    print(f"  {rlab(k)}  {dens:7.2f}   {np.median(dv):7.3f}  {np.median(np.ptp(avin[nb_e[sel]], axis=1)):7.3f}"
          f"   {np.percentile(rr, 10):5.2f}-{np.percentile(rr, 90):5.2f}      {ov:5.2f}     {sig[0]:.4f} / {sig[1]:.4f}  {sel.sum():5d}")

# ---- 4. sign ------------------------------------------------------------------------------------------
big = N.rec & (np.abs(N.ev) > 0.1)
print(f"\n[4] sign: recovered ASTs with |input - output| > 0.1 in F606W: {big.sum()}, of which "
      f"{np.mean(N.ev[big] > 0):.3f} came out brighter (input - output > 0); the model applies "
      f"observed = true {'+' if S.martin_sign else '-'} (input - output), so blends move "
      f"{'fainter (WRONG)' if S.martin_sign else 'brighter (correct)'}")

# ---- 5. zero points: AST inputs vs the observed ridge ---------------------------------------------------
av_app, ai_app = avin + zp["v"], a["InputF814W"].values.astype(float) + zp["i"]
rid = C.rcut.ridge
ms = rid.mag > 19.2
o = np.argsort(rid.mag[ms])
rc_ms = lambda m: np.interp(m, rid.mag[ms][o], rid.color[ms][o])
rg = (rid.mag > lo) & (rid.mag < 18.3)
og = np.argsort(rid.mag[rg])
rc_rg = lambda m: np.interp(m, rid.mag[rg][og], rid.color[rg][og])
print("\n[5] AST inputs (+ zero points) minus the observed ridge colour, median per F606W bin")
for a_, b_, f in ((17.0, 17.6, rc_rg), (17.6, 18.2, rc_rg), (19.3, 19.8, rc_ms), (19.8, 20.3, rc_ms), (20.3, 20.6, rc_ms)):
    sel = (av_app > a_) & (av_app < b_)
    print(f"  F606W {a_:.1f}-{b_:.1f}: {np.median((av_app - ai_app)[sel] - f(av_app[sel])):+.4f} mag ({sel.sum()} ASTs)")

# ---- 6/7. radial closure and where the model is narrower -------------------------------------------------
M = ctx(ast_errors="independent", radial="observed")       # Martin's error model
models = {"production (radii from fit stars)": C, "radii from loose-cut photometry": R, "Martin errors + radii": M}
cats = {}
for name, X in models.items():
    cats[name] = X.model(RUN, AGE)[1](np.random.RandomState(S.crn_seed), MU, EV, with_radius=True)
fr = {k: v for k, v in zip(["v", "i", "r"], [fit.v.values, fit.i.values, fit_r])}
win = (fr["v"] > lo) & (fr["v"] < hi)
print(f"\n[6] radial closure at run {RUN}, {AGE} Gyr, mu {MU}, E {EV}: median radius (px) and fraction within 200 px")
hdr = "  F606W        obs: med  <200"
for name in cats:
    hdr += f" | {name[:22]:>22s}: med  <200  KS"
print(hdr)
from scipy.stats import ks_2samp
mbins = np.arange(lo, hi + 1e-6, 0.5)
for a_, b_ in zip(mbins[:-1], mbins[1:]):
    so = win & (fr["v"] > a_) & (fr["v"] <= b_)
    line = f"  {a_:5.2f}-{b_:5.2f}   {np.median(fr['r'][so]):6.0f} {np.mean(fr['r'][so] < 200):5.3f}"
    for name, cat in cats.items():
        sc = (cat[:, 0] > a_) & (cat[:, 0] <= b_)
        line += (f" | {'':22s}  {np.median(cat[sc, 2]):6.0f} {np.mean(cat[sc, 2] < 200):5.3f} "
                 f"{ks_2samp(fr['r'][so], cat[sc, 2]).statistic:.3f}")
    print(line)


def widths(v, col, r):
    """robust colour width about the observed ridge per (magnitude region, radius bin)."""
    out = {}
    for reg, (a_, b_, f) in {"RGB 17.0-18.2": (17.0, 18.2, rc_rg), "MS 19.4-20.0": (19.4, 20.0, rc_ms),
                              "MS 20.0-20.6": (20.0, 20.6, rc_ms)}.items():
        for k in range(len(RB) - 1):
            sel = (v > a_) & (v < b_) & (rbin(r) == k)
            out[(reg, k)] = (rsig(col[sel] - f(v[sel])) if sel.sum() > 50 else np.nan, sel.sum())
    return out


wo = widths(fr["v"], fr["v"] - fr["i"], fr["r"])
wm = {name: widths(cat[:, 0], cat[:, 0] - cat[:, 1], cat[:, 2]) for name, cat in cats.items()}
print("\n[7] robust colour width (mag) about the observed ridge: observed vs each model, by radius")
for reg in ("RGB 17.0-18.2", "MS 19.4-20.0", "MS 20.0-20.6"):
    print(f"  {reg}")
    for k in range(len(RB) - 1):
        s_o, n_o = wo[(reg, k)]
        if not np.isfinite(s_o):
            continue
        line = f"    r {rlab(k)}: obs {s_o:.4f} ({n_o:5d})"
        for name in cats:
            s_m = wm[name][(reg, k)][0]
            line += f" | {name[:12]} {s_m:.4f} (x{s_o / s_m:.2f})"
        print(line)

# ---- 8. data only: is the extra width spatially coherent (differential reddening) or along y (CTE)? -----
sel = win & (fr["v"] > 19.4) & (fr["v"] < 20.6)
res = (fr["v"] - fr["i"])[sel] - rc_ms(fr["v"][sel])
fx, fy = fit.x.values[sel], fit.y.values[sel]
nb = 12
xe, ye = np.linspace(fx.min(), fx.max(), nb + 1), np.linspace(fy.min(), fy.max(), nb + 1)
med, err, cnt = np.full((nb, nb), np.nan), np.full((nb, nb), np.nan), np.zeros((nb, nb), int)
ix, iy = np.clip(np.digitize(fx, xe) - 1, 0, nb - 1), np.clip(np.digitize(fy, ye) - 1, 0, nb - 1)
for p in range(nb):
    for q in range(nb):
        m_ = (ix == p) & (iy == q)
        cnt[q, p] = m_.sum()
        if m_.sum() >= 40:
            med[q, p] = np.median(res[m_])
            err[q, p] = 1.25 * rsig(res[m_]) / np.sqrt(m_.sum())
ok = np.isfinite(med)
chi2 = np.nansum(((med - np.nanmedian(med)) / err) ** 2) / ok.sum()
cx, cy = (xe[:-1] + xe[1:]) / 2, (ye[:-1] + ye[1:]) / 2
CX, CY = np.meshgrid(cx, cy)
cr = np.hypot(CX - xc, CY - yc)
print(f"\n[8] upper-MS colour offset from the ridge in {nb}x{nb} cells ({ok.sum()} cells with >= 40 stars)")
print(f"  scatter of the cell medians {np.nanstd(med):.4f} mag vs their median noise {np.nanmedian(err):.4f}; "
      f"reduced chi2 about a constant {chi2:.1f}")
print(f"  correlation of the cell medians with y (CTE direction) {np.corrcoef(CY[ok], med[ok])[0, 1]:+.2f}, "
      f"with x {np.corrcoef(CX[ok], med[ok])[0, 1]:+.2f}, with radius {np.corrcoef(cr[ok], med[ok])[0, 1]:+.2f}")

# ---- figures ------------------------------------------------------------------------------------------
fig, axs = plt.subplots(2, 3, figsize=(17, 10), constrained_layout=True)
ax = axs[0, 0]
sub = rs.choice(len(inwin), min(25000, len(inwin)), replace=False)
ax.scatter((inwin.v - inwin.i).values[sub], inwin.v.values[sub], s=0.3, c=np.where(inwin.kept.values[sub], "0.6", "#d6604d"), lw=0)
ax.plot(rid.color, rid.mag, color="k", lw=0.8)
ax.set_xlim(0.3, 1.0); ax.set_ylim(hi + 0.1, lo - 0.1)
ax.set_title("loose-cut stars in the window: grey kept by the ridge cut, red cut", fontsize=10)
ax.set_xlabel("F606W - F814W"); ax.set_ylabel("F606W")
ax = axs[0, 1]
for name, cat in cats.items():
    sc = (cat[:, 0] > 19.4) & (cat[:, 0] < 20.6)
    ax.plot(np.sort(cat[sc, 2]), np.linspace(0, 1, sc.sum()), lw=1.2, label=name)
so = win & (fr["v"] > 19.4) & (fr["v"] < 20.6)
ax.plot(np.sort(fr["r"][so]), np.linspace(0, 1, so.sum()), color="k", lw=2, label="observed fit stars")
ax.set_xscale("log"); ax.set_xlim(20, 4000); ax.set_xlabel("radius (px)"); ax.set_ylabel("cumulative fraction")
ax.set_title("radial closure, F606W 19.4-20.6 (after every cut)", fontsize=10); ax.legend(fontsize=7)
ax = axs[0, 2]
sub = rs.choice(len(a), 20000, replace=False)
ax.scatter((av_app - ai_app)[sub], av_app[sub], s=0.5, color="#2166ac", lw=0, label="AST inputs + zero points")
ax.plot(rid.color, rid.mag, color="#d6604d", lw=1.2, label="observed ridge")
ax.set_xlim(0.3, 1.0); ax.set_ylim(hi + 0.5, lo - 0.5); ax.legend(fontsize=7)
ax.set_title("zero points: AST inputs on the observed ridge", fontsize=10); ax.set_xlabel("F606W - F814W")
ax = axs[1, 0]
if rows3:
    t = np.array(rows3)
    ax.plot(t[:, 0] + 1, t[:, 1], "o-", label="median |dF606W| to neighbours")
    ax.plot(t[:, 0] + 1, t[:, 2], "s-", label="F606W span of the 30")
    ax.plot(t[:, 0] + 1, t[:, 4] - t[:, 3], "^-", label="r_nb / r, p90 - p10")
    ax.set_xscale("log"); ax.set_xlabel("radius bin start (px)"); ax.legend(fontsize=7)
    ax.set_title("30-neighbour reach vs radius", fontsize=10)
ax = axs[1, 1]
for j, reg in enumerate(("RGB 17.0-18.2", "MS 19.4-20.0", "MS 20.0-20.6")):
    ks_ = [k for k in range(len(RB) - 1) if np.isfinite(wo[(reg, k)][0])]
    x = [(RB[k] + RB[k + 1]) / 2 for k in ks_]
    ax.plot(x, [wo[(reg, k)][0] / wm["production (radii from fit stars)"][(reg, k)][0] for k in ks_], "o-", label=reg)
ax.axhline(1, color="0.4", ls=":")
ax.set_xscale("log"); ax.set_xlabel("radius (px)"); ax.set_ylabel("observed / model colour width")
ax.set_title("where the data are wider than the production model", fontsize=10); ax.legend(fontsize=7)
ax = axs[1, 2]
lim = np.nanmax(np.abs(med - np.nanmedian(med)))
im = ax.pcolormesh(xe, ye, med - np.nanmedian(med), cmap="RdBu_r", vmin=-lim, vmax=lim)
ax.plot(xc, yc, "k+", ms=12)
fig.colorbar(im, ax=ax, label="median colour offset (mag), red = redder")
ax.set_aspect("equal"); ax.set_xlabel("x (px)"); ax.set_ylabel("y (px)")
ax.set_title(f"upper-MS colour offset map: reduced chi2 {chi2:.1f}", fontsize=10)
fig.suptitle(f"M15 sCMD error-model diagnostics (run {RUN}, {AGE} Gyr, mu {MU}, E(V-I) {EV})")
fig.savefig(f"{OUT}/diagnostics.png", dpi=110)
print(f"\n{OUT}/diagnostics.png")
