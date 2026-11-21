"""Cheap checks of the context: counts, ridge settings, centres, AST list, ties, seeds."""
import sys, json, time, numpy as np, pandas as pd
sys.path.insert(0, '/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox')
from isochronetoolbox.scmd import Setup, context
from isochronetoolbox.scmd.pipeline import _window, seed, select_age
from isochronetoolbox.scmd.ast import centre
G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
S = Setup(**{k: v for k, v in json.load(open(G + '/setup.json')).items() if k in Setup.__dataclass_fields__})
print("ks_percentile (not in setup.json -> default):", S.ks_percentile, " ks_3sigma:", S.ks_3sigma)
t = time.time(); ctx = context(S); print(f"context {time.time()-t:.1f}s")
rc = ctx.rcut
print("n_target", ctx.n_target, " len(obs)", len(ctx.obs), " len(obs3)", len(ctx.obs3), " v_window", ctx.v_window)
print("rcut: n_trace", rc.n_trace, " scales", rc.scales, " bw", rc.bw, " width", rc.width, " margin", rc.margin,
      " mag_range", rc.mag_range, " color_range", rc.color_range, " exclude", rc.exclude)
# raw file: VIvega vs Vvega - Ivega; rounding
raw = pd.read_csv(S.obs_raw, sep=r"\s+")
print("raw columns", list(raw.columns)[:20], " n_raw", len(raw))
dvi = raw["VIvega"] - (raw["Vvega"] - raw["Ivega"])
print("VIvega - (Vvega-Ivega): max|.|", float(np.abs(dvi).max()), " frac nonzero(>1e-6)", float((np.abs(dvi) > 1e-6).mean()))
o3 = ctx.obs3
print("obs3 distinct colours", len(np.unique(o3[:, 0])), " distinct mags", len(np.unique(o3[:, 1])),
      " distinct (c,m) pairs", len(np.unique(o3, axis=0)))
# centres: fit stars vs ASTs
a = pd.read_csv(S.ast, dtype=np.float32)
print("AST columns", list(a.columns), " n_ast", len(a))
xc_a, yc_a = centre(a["InputX"].values, a["InputY"].values)
fit = raw.loc[:, ["x", "y"]]
# fit stars used for the radial distribution = ridge-selected observed stars
obs_sel = pd.DataFrame({"x": raw["x"], "y": raw["y"]})
print("AST centre", (round(xc_a, 2), round(yc_a, 2)))
# reconstruct the fit-star centre the way ArtificialStars does (from ctx.ast.dist we only have radii); recompute:
from isochronetoolbox.scmd.ridge import RidgeCut
v, vi = raw["Vvega"].values, raw["VIvega"].values
keep = rc.select(vi, v, rc.ridge)
xs, ys = raw["x"].values[keep], raw["y"].values[keep]
xc_f, yc_f = centre(xs, ys)
print("fit-star centre", (round(xc_f, 2), round(yc_f, 2)), " offset AST-fit (px):", round(float(np.hypot(xc_a-xc_f, yc_a-yc_f)), 2),
      " n fit stars", keep.sum())
print("radial distribution: n", len(ctx.ast.dist), " min/median/max r", ctx.ast.dist.min(), np.median(ctx.ast.dist), ctx.ast.dist.max())
# AST input radius range
r_ast = np.hypot(a["InputX"].values - xc_a, a["InputY"].values - yc_a)
print("AST radius min/median/max", r_ast.min(), np.median(r_ast), r_ast.max())
ne = ctx.ast.nearest
print("AST recovered fraction", ne.rec.mean())
# input magnitudes in VEGA
vin = a["InputF606W"].values + S.zp_v; iin = a["InputF814W"].values + S.zp_i
print("AST input V range (VEGA)", vin.min(), vin.max(), " input V-I range", (vin-iin).min(), (vin-iin).max())
inwin = (vin > 15.02) & (vin < 19.02)
print("ASTs with input V in window:", inwin.sum())
# halves
ha, hb = ctx.ast.halves()
for nm, h in (("A", ha), ("B", hb)):
    aa = pd.read_csv(S.ast, dtype=np.float32)
print("half sizes", ha.n_ast, hb.n_ast)
half = np.random.RandomState(0).random_sample(len(a)) < 0.5
for nm, m in (("A", half), ("B", ~half)):
    xc, yc = centre(a["InputX"].values[m], a["InputY"].values[m])
    print(f"half {nm} AST centre", (round(xc, 2), round(yc, 2)), " shift vs full", round(float(np.hypot(xc-xc_a, yc-yc_a)), 2))
print("radial dist identical in halves:", np.array_equal(ha.dist, ctx.ast.dist), np.array_equal(hb.dist, ctx.ast.dist))
# grid seed and the cell's position in the stream
sd = seed(10099, 11.4)
print("grid seed(10099, 11.4) =", sd, " (< 10000 would collide with a realisation seed)")
dm = S.dm; ev = S.ebv
idx = [(i, j) for i, d in enumerate(dm) for j, e in enumerate(ev) if abs(d-13.38) < 5e-4 and abs(e-0.015) < 5e-4]
print("dm grid n", len(dm), " evi grid n", len(ev), " cell index (dm, evi)", idx, " -> position in stream", idx[0][0]*len(ev)+idx[0][1])
print("dm grid values near 13.38:", dm[20:26], " evi grid:", ev)
# ties in masses
iso = ctx.isochrones(10099)
age_yr, rows = select_age(iso, 11.4, S.exact_age)
m = rows["initial_mass"].values
print("isochrone age", age_yr, " points", len(m), " unique masses", len(np.unique(m)), " monotonic increasing?", bool(np.all(np.diff(m) > 0)))
