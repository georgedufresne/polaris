"""Off-grid penalty at the best cell: KS vs data at sub-step offsets, common random numbers (one seed)."""
import sys, json, numpy as np
sys.path.insert(0, '/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox')
from isochronetoolbox.scmd import Setup, context
from isochronetoolbox.scmd.pipeline import _window
from isochronetoolbox.scmd.stats import ks_2d
axis = sys.argv[1]
G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
S = Setup(**{k: v for k, v in json.load(open(G + '/setup.json')).items() if k in Setup.__dataclass_fields__})
ctx = context(S)
_, make = ctx.model(10099, 11.4)
if axis == "evi":
    pts = [(13.38, e) for e in (0.0125, 0.0140, 0.0150, 0.0160, 0.0175)]
else:
    pts = [(m, 0.015) for m in (13.370, 13.375, 13.385, 13.390, 13.400)]
for sd in (11, 12):
    for dm, e in pts:
        p = ctx.points(make(np.random.RandomState(sd), dm, e)); p = p[_window(p, S.ks_window)]
        print(f"seed {sd} mu {dm:.3f} E(V-I) {e:.4f}: KS {ks_2d(ctx.obs3, p, ctx.obs3_cdf)[0]:.5f}", flush=True)
