"""Shared helpers for the ast-correction audit (read-only use of the toolbox)."""
import sys, json, time, types
import numpy as np, pandas as pd
sys.path.insert(0, '/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox')
from isochronetoolbox.scmd import Setup, context
from isochronetoolbox.scmd.pipeline import _window, regenerate
from isochronetoolbox.scmd.stats import ks_2d, rank2d
from isochronetoolbox.scmd.ast import ArtificialStars, NearestErrors, centre
from scipy.interpolate import LinearNDInterpolator

G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
RUN, AGE, DM, EVI = 10099, 11.4, 13.38, 0.015
OUT = '/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/ast-correction'


def setup():
    return Setup(**{k: v for k, v in json.load(open(G + '/setup.json')).items() if k in Setup.__dataclass_fields__})


def window_pts(ctx, cat):
    pts = ctx.points(cat)
    return pts[_window(pts, ctx.s.ks_window)]


def sample(ctx, make, ast, seed):
    """One sCMD of the best cell drawn with AST model ``ast`` and RandomState(seed), cut to the KS window."""
    full = ctx.ast
    try:
        ctx.ast = ast
        cat = make(np.random.RandomState(int(seed)), DM, EVI)
    finally:
        ctx.ast = full
    return window_pts(ctx, cat)


class Ref:
    def __init__(self, pts, name):
        self.pts, self.name = pts, name
        self.cdf = rank2d(pts) / len(pts)


def ks_multi(sim, refs, q=95.0):
    """KS (as stats.ks_2d) of one simulated set against several references, building
    the simulated ECDF and its interpolator once."""
    sim_cdf = rank2d(sim) / len(sim)
    f = LinearNDInterpolator(sim, sim_cdf, fill_value=0)
    out = []
    for r in refs:
        d = np.abs(r.cdf - f(r.pts))
        out.append(float(np.percentile(d, q)))
    return out
