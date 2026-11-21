#!/usr/bin/env python3
"""Error-model and photometry diagnostics for one cluster (the M15 checks, generalised).

  python diag_generic.py CONFIG.json

CONFIG: {"dir": experiment dir (code/, setup.json), "label": str, "cell": [run, age, mu, E],
         "overrides": {setup overrides for the model}, "ms": [F606W lo, hi], "rgb": [lo, hi],
         "sgb_color": [lo, hi], "sgb_mag": [lo, hi]}
Prints: [1] ridge-cut efficiency by radius, [2] completeness floor, [3] colour width obs vs model by
radius, [4] spatial colour-offset map, [5] per-star SGB-vs-local-colour slope -> F606W share k of the
spatial pattern (reddening k ~ 2.8, F606W zero point 1, F814W zero point 0)."""
import json
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CFG = json.load(open(sys.argv[1]))
D = CFG["dir"]
sys.path.insert(0, f"{D}/code")
from isochronetoolbox.scmd import Setup, context                # noqa: E402

RB = np.array([0, 100, 200, 400, 700, 1200, 2000, 4500])
base = {k: v for k, v in json.load(open(f"{D}/setup.json")).items() if k in Setup.__dataclass_fields__}


def ctx(**kw):
    d = dict(base); d.update(CFG.get("overrides", {})); d.update(kw)
    return context(Setup(**d))


def rsig(x):
    x = x[np.isfinite(x)]
    return 1.4826 * np.median(np.abs(x - np.median(x))) if len(x) > 50 else np.nan


C = ctx()
R = ctx(radial_source="raw")
S = C.s
lo, hi = C.v_window
fit = C.ast._args[0]
xc, yc = C.ast.centre_xy
v, c, x, y = fit.v.values, (fit.v - fit.i).values, fit.x.values, fit.y.values
r = np.hypot(x - xc, y - yc)
rid = C.rcut.ridge
ms_lo, ms_hi = CFG["ms"]
m1 = rid.mag > ms_lo - 0.2
o1 = np.argsort(rid.mag[m1])
ms_c = lambda V: np.interp(V, rid.mag[m1][o1], rid.color[m1][o1])
g_lo, g_hi = CFG["rgb"]
m3 = (rid.mag > g_lo - 0.2) & (rid.mag < g_hi + 0.2)
o3 = np.argsort(rid.mag[m3])
rg_c = lambda V: np.interp(V, rid.mag[m3][o3], rid.color[m3][o3])
sc_lo, sc_hi = CFG["sgb_color"]
sm_lo, sm_hi = CFG["sgb_mag"]
m2 = (rid.mag > sm_lo) & (rid.mag < sm_hi) & (rid.color > sc_lo - 0.02) & (rid.color < sc_hi + 0.02)
o2 = np.argsort(rid.color[m2])
sgb_v = lambda col: np.interp(col, rid.color[m2][o2], rid.mag[m2][o2])
print(f"{CFG['label']}: window {lo:.2f}-{hi:.2f}, fit stars {len(fit)}, centre ({xc:.0f}, {yc:.0f})")

# [1] ridge-cut efficiency in the core
loose = R.radial_stars.copy()
inwin = loose[(loose.v > lo) & (loose.v < hi)].reset_index(drop=True)
inwin = inwin.merge(fit[["x", "y"]].drop_duplicates().assign(kept=True), on=["x", "y"], how="left")
inwin["kept"] = inwin["kept"].fillna(False).astype(bool)
ir = np.hypot(inwin.x - xc, inwin.y - yc).values
print("\n[1] ridge-cut efficiency (kept / loose-cut stars) by radius: " + ", ".join(
    f"{RB[k]}-{RB[k + 1]} px {inwin.kept[(ir >= RB[k]) & (ir < RB[k + 1])].mean():.3f}"
    for k in range(len(RB) - 1) if ((ir >= RB[k]) & (ir < RB[k + 1])).sum() > 30))

# [2] completeness floor
cc = C.ast.nearest.completeness(v - S.zp_v, (v - c) - S.zp_i, r)
print(f"[2] fit stars with completeness < 0.1: {np.mean(cc < 0.1):.4f}; < 0.2: {np.mean(cc < 0.2):.4f}; "
      f"median in the inner 100 px {np.median(cc[r < 100]) if (r < 100).any() else np.nan:.3f}")

# [3] widths by radius, observed vs model
run, age, mu, ev = CFG["cell"]
cats = {name: X.model(run, age)[1](np.random.RandomState(S.crn_seed), mu, ev, with_radius=True)
        for name, X in (("model", C), ("model, raw radii", R))}
print(f"\n[3] robust colour width about the observed ridge (model: run {run}, {age} Gyr, mu {mu}, E {ev})")
rows3 = {}
for reg, (a_, b_, f) in {f"RGB {g_lo}-{g_hi}": (g_lo, g_hi, rg_c), f"MS {ms_lo}-{ms_hi}": (ms_lo, ms_hi, ms_c)}.items():
    print(f"  {reg}:   r (px)      obs    model  ratio  missing sigma")
    for k in range(len(RB) - 1):
        so = (v > a_) & (v < b_) & (r >= RB[k]) & (r < RB[k + 1])
        cat = cats["model"]
        sm = (cat[:, 0] > a_) & (cat[:, 0] < b_) & (cat[:, 2] >= RB[k]) & (cat[:, 2] < RB[k + 1])
        wo, wm = rsig(c[so] - f(v[so])), rsig((cat[sm, 0] - cat[sm, 1]) - f(cat[sm, 0]))
        if np.isfinite(wo) and np.isfinite(wm):
            rows3[(reg, k)] = (wo, wm)
            print(f"           {RB[k]:5d}-{RB[k + 1]:5d}  {wo:.4f}  {wm:.4f}  {wo / wm:5.2f}   {np.sqrt(max(wo * wo - wm * wm, 0)):.4f}")

# [4] spatial colour-offset map (upper MS, r > 0)
sel = (v > ms_lo) & (v < ms_hi)
res = c[sel] - ms_c(v[sel])
nb = 12
fx, fy = x[sel], y[sel]
xe, ye = np.linspace(fx.min(), fx.max(), nb + 1), np.linspace(fy.min(), fy.max(), nb + 1)
med, err = np.full((nb, nb), np.nan), np.full((nb, nb), np.nan)
ix, iy = np.clip(np.digitize(fx, xe) - 1, 0, nb - 1), np.clip(np.digitize(fy, ye) - 1, 0, nb - 1)
for p in range(nb):
    for q in range(nb):
        mm = (ix == p) & (iy == q)
        if mm.sum() >= 40:
            med[q, p] = np.median(res[mm]); err[q, p] = 1.25 * rsig(res[mm]) / np.sqrt(mm.sum())
ok = np.isfinite(med)
chi2 = np.nansum(((med - np.nanmedian(med)) / err) ** 2) / ok.sum()
print(f"\n[4] upper-MS colour offset map ({ok.sum()} cells): cell-median scatter {np.nanstd(med):.4f} vs noise "
      f"{np.nanmedian(err):.4f}; reduced chi2 {chi2:.1f}")

# [5] per-star test: SGB F606W offset vs local MS colour offset
h = 0.01
MS = sel & (r > 300)
SG = (c > sc_lo) & (c < sc_hi) & (v > sm_lo) & (v < sm_hi) & (np.abs(v - sgb_v(c)) < 0.2) & (r > 400)
s_ms = np.median((ms_c(v[MS] + h) - ms_c(v[MS] - h)) / (2 * h))
s_sg = np.median((sgb_v(c[SG] + h) - sgb_v(c[SG] - h)) / (2 * h))
ms_res = c[MS] - ms_c(v[MS])
tree = cKDTree(np.column_stack([x[MS], y[MS]]))
print(f"\n[5] {SG.sum()} SGB stars, {MS.sum()} MS stars; ridge slopes dcol/dV (MS) {s_ms:.4f}, dV/dcol (SGB) {s_sg:.2f}")
pred = {k: (k - s_sg) / (1 - s_ms * k) for k in (2.8, 1.0, 0.0)}
print("    predicted slope: reddening {:.2f}, F606W zero point {:.2f}, F814W zero point {:.2f}".format(pred[2.8], pred[1.0], pred[0.0]))
rs = np.random.RandomState(0)
for K in (100, 150, 300):
    _, idx = tree.query(np.column_stack([x[SG], y[SG]]), k=K)
    dc = np.median(ms_res[idx], axis=1)
    dv = v[SG] - sgb_v(c[SG])

    def fitline(ii):
        q = np.quantile(dc[ii], np.linspace(0, 1, 9))
        b = np.clip(np.digitize(dc[ii], q[1:-1]), 0, 7)
        X = np.array([np.median(dc[ii][b == j]) for j in range(8)])
        Y = np.array([np.median(dv[ii][b == j]) for j in range(8)])
        return np.polyfit(X, Y, 1)[0]
    sl = fitline(np.arange(len(dc)))
    bs = np.array([fitline(rs.randint(0, len(dc), len(dc))) for _ in range(500)])
    noise = 1.25 * rsig(ms_res) / np.sqrt(K)
    att = (np.var(dc) - noise ** 2) / np.var(dc)
    slc = sl / att
    k_est = (slc + s_sg) / (1 + s_ms * slc)
    print(f"    K {K}: local colour offset rms {np.std(dc):.4f} (noise {noise:.4f}); slope {sl:.2f} +/- {bs.std():.2f}, "
          f"noise-corrected {slc:.2f} -> F606W share k = {k_est:.2f}")

fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), constrained_layout=True)
ax = axs[0]
for reg in sorted({k[0] for k in rows3}):
    ks_ = [k for k in range(len(RB) - 1) if (reg, k) in rows3]
    ax.plot([(RB[k] + RB[k + 1]) / 2 for k in ks_], [rows3[(reg, k)][0] / rows3[(reg, k)][1] for k in ks_], "o-", label=reg)
ax.axhline(1, color="0.4", ls=":"); ax.set_xscale("log"); ax.legend(fontsize=8)
ax.set_xlabel("radius (px)"); ax.set_ylabel("observed / model colour width")
ax.set_title("where the data are wider than the model", fontsize=10)
ax = axs[1]
lim = np.nanmax(np.abs(med - np.nanmedian(med)))
im = ax.pcolormesh(xe, ye, med - np.nanmedian(med), cmap="RdBu_r", vmin=-lim, vmax=lim)
ax.plot(xc, yc, "k+", ms=12); fig.colorbar(im, ax=ax, label="median colour offset (mag), red = redder")
ax.set_aspect("equal"); ax.set_title(f"upper-MS colour offset map: reduced chi2 {chi2:.1f}", fontsize=10)
fig.suptitle(f"{CFG['label']}: photometry and error-model diagnostics")
out = f"{D}/diag/diagnostics.png"
fig.savefig(out, dpi=110); print(out)
