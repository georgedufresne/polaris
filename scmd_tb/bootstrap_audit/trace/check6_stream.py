"""(1) Does the random-stream consumption of one sCMD depend on the AST list or the ridge?  If not, the
regenerated reference (cell 256 of the (run, age) stream) has the same masses, binaries and radii as the
grid's sCMD, whatever list it takes its errors from.  (2) Extent of the k=30 AST neighbourhoods."""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, '/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox')
from isochronetoolbox.scmd import Setup, context
from isochronetoolbox.scmd.pipeline import seed
G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
S = Setup(**{k: v for k, v in json.load(open(G + '/setup.json')).items() if k in Setup.__dataclass_fields__})
ctx = context(S)
full = ctx.ast; A, B = full.halves()
_, make = ctx.model(10099, 11.4)
calls = {"select": 0, "trace_box": 0}
o_sel, o_tb = ctx.rcut.select, ctx.rcut.in_trace_box
ctx.rcut.select = lambda *a, **k: (calls.__setitem__("select", calls["select"] + 1), o_sel(*a, **k))[1]
ctx.rcut.in_trace_box = lambda *a, **k: (calls.__setitem__("trace_box", calls["trace_box"] + 1), o_tb(*a, **k))[1]
cells = [(13.15, 0.0), (13.15, 0.005), (13.38, 0.015), (13.45, 0.05)]
for nm, ast in (("full", full), ("half A", A), ("half B", B)):
    ctx.ast = ast
    rs = np.random.RandomState(seed(10099, 11.4))
    calls.update(select=0, trace_box=0)
    cats = [make(rs, d, e) for d, e in cells]
    nxt = rs.random_sample(3)
    print(f"{nm:7s}: next uniforms after {len(cells)} cells {np.round(nxt, 12)}  pool batches {calls['trace_box']}  select calls {calls['select']}"
          f"  (1 select per sCMD = no top-up draws)")
ctx.ast = full
ctx.rcut.select, ctx.rcut.in_trace_box = o_sel, o_tb
# (2) AST neighbourhood extent
ne = full.nearest
a = pd.read_csv(S.ast, dtype=np.float32)
from isochronetoolbox.scmd.ast import centre
xc, yc = centre(a["InputX"].values, a["InputY"].values)
r = np.hypot(a["InputX"].values - xc, a["InputY"].values - yc)
V = a["InputF606W"].values + S.zp_v
win = (V > 15.02) & (V < 19.02)
for lo, hi in ((0, 150), (150, 300), (300, 600), (600, 1200), (1200, 3000)):
    m = win & (r >= lo) & (r < hi)
    nb = ne.nbr[m]
    dr = np.abs(r[nb] - r[m][:, None]); dv = np.abs(V[nb] - V[m][:, None])
    print(f"r {lo:4d}-{hi:4d} px: ASTs in window {m.sum():5d}  neighbours' |dr| median {np.median(dr):6.1f} p90 {np.percentile(dr, 90):6.1f} px;"
          f" |dV| median {np.median(dv):.3f} p90 {np.percentile(dv, 90):.3f} mag; recovered {ne.rec[m].mean():.3f}")
# where do synthetic radii fall? fraction of fit stars at r<150, 150-300
d = full.dist
print("observed fit-star radii: frac r<150 %.4f, 150-300 %.4f, 300-600 %.4f" % ((d < 150).mean(), ((d >= 150) & (d < 300)).mean(), ((d >= 300) & (d < 600)).mean()))
print("AST radii (all):         frac r<150 %.4f, 150-300 %.4f, 300-600 %.4f" % ((r < 150).mean(), ((r >= 150) & (r < 300)).mean(), ((r >= 300) & (r < 600)).mean()))
