"""Cheap checks: rank2d tie handling, convex-hull coverage, pool size for tracing, AST criteria."""
import sys, json, time, numpy as np, pandas as pd
sys.path.insert(0, '/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox')
from isochronetoolbox.scmd import Setup, context
from isochronetoolbox.scmd.pipeline import _window
from isochronetoolbox.scmd.stats import ks_2d, rank2d
from scipy.spatial import Delaunay
G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
S = Setup(**{k: v for k, v in json.load(open(G + '/setup.json')).items() if k in Setup.__dataclass_fields__})
ctx = context(S)
o3 = ctx.obs3; N = len(o3)

# 1. rank2d vs brute force on the observed (tied) points
rs = np.random.RandomState(1)
sub = rs.choice(N, 400, replace=False)
r2 = rank2d(o3)
strict = np.array([np.sum((o3[:, 0] < o3[k, 0]) & (o3[:, 1] < o3[k, 1])) for k in sub])
le = np.array([np.sum((o3[:, 0] <= o3[k, 0]) & (o3[:, 1] <= o3[k, 1])) - 1 for k in sub])   # minus self
d_strict = (r2[sub] - strict) / N
print("rank2d on observed: (rank2d - strict)/N  mean %.5f  p95 %.5f  max %.5f ; (le - strict)/N p95 %.5f max %.5f" % (
    d_strict.mean(), np.percentile(np.abs(d_strict), 95), np.abs(d_strict).max(),
    np.percentile((le - strict) / N, 95), ((le - strict) / N).max()))
print("  rank2d within [strict, le]:", bool(np.all((r2[sub] >= strict) & (r2[sub] <= le))))
# tiny tie example
pts = np.array([[0.5, 17.0], [0.5, 16.0], [0.5, 15.0], [0.4, 15.5]])
print("  toy: points", pts.tolist(), " rank2d", rank2d(pts).tolist(), " strict", [int(np.sum((pts[:, 0] < p[0]) & (pts[:, 1] < p[1]))) for p in pts])

# 2. pool size used for tracing (record what rcut.trace receives)
seen = []
orig = ctx.rcut.trace
def rec(c, m):
    seen.append(len(c)); return orig(c, m)
ctx.rcut.trace = rec
_, make = ctx.model(10099, 11.4)
t = time.time()
cats = [make(np.random.RandomState(i), 13.38, 0.015) for i in (101, 102)]
print("two sCMDs %.1fs; trace input sizes (latest code: pool[:n_trace])" % (time.time() - t), seen, " n_trace", ctx.rcut.n_trace)
ctx.rcut.trace = orig
# how big is the pool itself? replicate fill's first batch rate with the same model: count pool by monkeypatching in_trace_box
cnt = []
orig_itb = ctx.rcut.in_trace_box
def itb(c, m):
    k = orig_itb(c, m); cnt.append(int(k.sum())); return k
ctx.rcut.in_trace_box = itb
_ = make(np.random.RandomState(103), 13.38, 0.015)
ctx.rcut.in_trace_box = orig_itb
print("  in-trace-box stars per fill batch (before recovery mask is ANDed):", cnt)

# 3. convex-hull coverage: observed points outside a model sCMD's hull vs model reference outside another model's hull
P = [ctx.points(c) for c in cats]
P = [p[_window(p, S.ks_window)] for p in P]
for nm, ref, sim in (("data vs model", o3, P[0]), ("model ref vs model", P[1], P[0])):
    tri = Delaunay(sim)
    out = tri.find_simplex(ref) < 0
    ref_cdf = rank2d(ref) / len(ref)
    ks, d = ks_2d(ref, sim, ref_cdf)
    print(f"{nm}: n_ref {len(ref)} n_sim {len(sim)}  outside hull {out.sum()} ({100*out.mean():.3f}%)  "
          f"KS(p95) {ks:.5f}  max|d| {d.max():.4f}  mean d outside hull {d[out].mean() if out.any() else 0:.4f}  "
          f"p95 without hull-outside points {np.percentile(d[~out], 95):.5f}")
np.save('/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/trace/two_scmds.npy', np.array(P, dtype=object), allow_pickle=True)

# 4. AST recovery criteria: found-but-failed fraction
a = pd.read_csv(S.ast, dtype=np.float32)
vin, iin, vout, iout = (a[c].values.astype(float) for c in ("InputF606W", "InputF814W", "OutputF606W", "OutputF814W"))
dpos = np.hypot(a["InputX"].values - a["OutputX"].values, a["InputY"].values - a["OutputY"].values)
print("AST output sentinel check: OutputF606W min/max", vout.min(), vout.max(), " unique top values", np.unique(np.round(vout, 3))[-3:])
found = np.isfinite(vout) & (np.abs(vout) < 50) & (np.abs(vout - vin) < 5)
rec = (dpos < 0.5) & (np.abs(vin - vout) < 0.75) & (np.abs(iin - iout) < 0.75)
win = (vin + S.zp_v > 15.02) & (vin + S.zp_v < 19.02)
print("in window: n", win.sum(), " recovered", (rec & win).sum(), f"({(rec & win).sum()/win.sum():.3f})",
      " found (|out-in|<5) but not recovered", (found & ~rec & win).sum(), f"({(found & ~rec & win).sum()/win.sum():.4f})")
fr = found & ~rec & win
if fr.any():
    print("  of those: dpos>=0.5 & |dV|,|dI|<0.75:", int((fr & (dpos >= 0.5) & (np.abs(vin-vout) < 0.75) & (np.abs(iin-iout) < 0.75)).sum()),
          "  |dV| or |dI| >= 0.75:", int((fr & ((np.abs(vin-vout) >= 0.75) | (np.abs(iin-iout) >= 0.75))).sum()),
          "  median dpos", float(np.median(dpos[fr])), " median (in-out)V", float(np.median((vin-vout)[fr])))
# AST input colour coverage vs V
V = vin + S.zp_v; VI = (vin + S.zp_v) - (iin + S.zp_i)
for lo, hi in ((15.0, 16.0), (16.0, 17.0), (17.0, 17.5), (17.5, 18.0), (18.0, 19.0)):
    m = (V >= lo) & (V < hi)
    om = (o3[:, 1] >= lo) & (o3[:, 1] < hi)
    print(f"V {lo}-{hi}: AST n {m.sum():5d} input V-I p2/p50/p98 {np.percentile(VI[m], 2):.3f}/{np.median(VI[m]):.3f}/{np.percentile(VI[m], 98):.3f} | "
          f"observed n {om.sum():5d} V-I p2/p50/p98 {np.percentile(o3[om, 0], 2):.3f}/{np.median(o3[om, 0]):.3f}/{np.percentile(o3[om, 0], 98):.3f}")
