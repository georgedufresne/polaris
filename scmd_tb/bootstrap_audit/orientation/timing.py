import sys, json, time, numpy as np
sys.path.insert(0, '/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox')
from isochronetoolbox.scmd import Setup, context
from isochronetoolbox.scmd.pipeline import _window
from isochronetoolbox.scmd.stats import ks_2d, rank2d
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import Delaunay
G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
t = time.time()
S = Setup(**{k: v for k, v in json.load(open(G + '/setup.json')).items() if k in Setup.__dataclass_fields__})
ctx = context(S)
print('context', time.time() - t, len(ctx.obs3)); sys.stdout.flush()
t = time.time()
_, make = ctx.model(10099, 11.4)
print('model', time.time() - t)
t = time.time()
pts = ctx.points(make(np.random.RandomState(1000), 13.38, 0.015)); pts = pts[_window(pts, S.ks_window)]
print('scmd', time.time() - t, len(pts))
t = time.time(); ks = ks_2d(ctx.obs3, pts, ctx.obs3_cdf)[0]; print('ks_2d', time.time() - t, ks)
t = time.time(); r = rank2d(pts); print('rank2d sim', time.time() - t)
t = time.time(); f = LinearNDInterpolator(pts, r / len(pts), fill_value=0)(ctx.obs3); print('LNDI', time.time() - t)
t = time.time(); tri = Delaunay(pts); print('Delaunay', time.time() - t)
t = time.time(); s = tri.find_simplex(ctx.obs3); print('find_simplex', time.time() - t, (s < 0).sum())
both = np.vstack([pts, ctx.obs3]); mask = np.r_[np.ones(len(pts), bool), np.zeros(len(ctx.obs3), bool)]
t = time.time(); rc = rank2d(both, mask)[len(pts):]; print('cross rank2d', time.time() - t)
print('ties: distinct colours', len(np.unique(ctx.obs3[:, 0])), 'distinct mags', len(np.unique(ctx.obs3[:, 1])))
