"""ECDF orientation: the 95th-percentile statistic in each of the four quadrant conventions (exact cross-counts),
data vs model and model vs model; and where in the CMD the top 5 % of |d| sit (toolbox convention)."""
import sys, json, numpy as np
sys.path.insert(0, '/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox')
from isochronetoolbox.scmd import Setup, context
from isochronetoolbox.scmd.pipeline import _window
from isochronetoolbox.scmd.stats import rank2d
G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
S = Setup(**{k: v for k, v in json.load(open(G + '/setup.json')).items() if k in Setup.__dataclass_fields__})
ctx = context(S)
full = ctx.ast; A, B = full.halves()
_, make = ctx.model(10099, 11.4)
def scmd(ast, sd):
    ctx.ast = ast
    try:
        p = ctx.points(make(np.random.RandomState(sd), 13.38, 0.015))
    finally:
        ctx.ast = full
    return p[_window(p, S.ks_window)]
def stat(ref, sim, sx, sy):
    f = np.array([sx, sy])
    r, s = ref * f, sim * f
    rc = rank2d(r) / len(r)
    pts = np.vstack([s, r]); mask = np.r_[np.ones(len(s), bool), np.zeros(len(r), bool)]
    ex = rank2d(pts, mask)[len(s):] / len(s)
    d = rc - ex
    return np.percentile(np.abs(d), 95), d, rc
m1, m2 = scmd(full, 300), scmd(full, 400)
ma, mb = scmd(A, 401), scmd(B, 301)
names = {(1, 1): "bluer & brighter (toolbox)", (1, -1): "bluer & fainter (M. Ying: V flipped)",
         (-1, 1): "redder & brighter", (-1, -1): "redder & fainter"}
for (sx, sy), nm in names.items():
    kd = stat(ctx.obs3, m1, sx, sy)[0]; kn = stat(m2, m1, sx, sy)[0]; ks_ = stat(ma, mb, sx, sy)[0]
    z = stat(ctx.obs3, m1, sx, sy)[2]
    print(f"{nm:38s}: data-vs-model {kd:.5f}  model-vs-model shared {kn:.5f}  split {ks_:.5f}  ratio data/shared {kd/kn:5.2f}"
          f"  obs points with ECDF=0: {np.mean(z == 0):.3f}")
k, d, rc = stat(ctx.obs3, m1, 1, 1)
top = np.abs(d) >= k
v = ctx.obs3[:, 1]; c = ctx.obs3[:, 0]
print("toolbox convention, data vs model: KS %.5f; top-5%% points: sign(data-model) + frac %.2f" % (k, np.mean(d[top] > 0)))
for lo, hi in ((15.02, 16.5), (16.5, 17.0), (17.0, 17.5), (17.5, 18.0), (18.0, 18.5), (18.5, 19.02)):
    m = (v >= lo) & (v < hi)
    print(f"  V {lo:5.2f}-{hi:5.2f}: stars {m.sum():6d}  share of top-5% {np.mean(top & m)/np.mean(top):.3f}  median ECDF {np.median(rc[m]):.3f}  median |d| {np.median(np.abs(d[m])):.4f}")
