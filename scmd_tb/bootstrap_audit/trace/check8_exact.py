"""Interpolated vs exact cross-count ECDF for the 95th-percentile statistic: data vs model and model vs model."""
import sys, json, numpy as np
sys.path.insert(0, '/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox')
from isochronetoolbox.scmd import Setup, context
from isochronetoolbox.scmd.pipeline import _window
from isochronetoolbox.scmd.stats import ks_2d, rank2d
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
def both(ref, sim, ref_cdf=None):
    ref_cdf = rank2d(ref) / len(ref) if ref_cdf is None else ref_cdf
    ks_i, d_i = ks_2d(ref, sim, ref_cdf)
    pts = np.vstack([sim, ref]); mask = np.r_[np.ones(len(sim), bool), np.zeros(len(ref), bool)]
    ex = rank2d(pts, mask)[len(sim):] / len(sim)
    d_e = np.abs(ref_cdf - ex)
    return ks_i, np.percentile(d_e, 95), np.percentile(np.abs(d_i - d_e), 95)
for k in range(3):
    m = scmd(full, 300 + k)
    print("data vs model(full) seed %d: p95 interp %.5f exact %.5f | p95 |interp-exact| %.5f" % ((300 + k,) + both(ctx.obs3, m, ctx.obs3_cdf)), flush=True)
    r = scmd(full, 400 + k)
    print("  model(full) vs model(full):  p95 interp %.5f exact %.5f | %.5f" % both(r, m), flush=True)
    ra, mb = scmd(A, 400 + k), scmd(B, 300 + k)
    print("  model(A) vs model(B):        p95 interp %.5f exact %.5f | %.5f" % both(ra, mb), flush=True)
