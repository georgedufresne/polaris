"""Shared helpers for the reference-selection audit (read-only use of the toolbox).

KS evaluation is the toolbox's ks_2d, restructured only for speed:
* the realisation's Delaunay triangulation is built once and evaluated at several
  reference point sets;
* the query (reference) points are sorted (lexsort by magnitude, then colour) so
  scipy's simplex walk starts near each point -- bit-identical interpolated values
  (checked against ks_2d in check_equivalence below), ~5x faster.
"""
import json
import os
import sys
import time

import numpy as np

W = "/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/reference-selection"
OUT = os.path.join(W, "out")
G = "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02"
LATEST = "/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox"
GRIDCODE = os.path.join(W, "gridcode")
RUN, AGE, MU, EVI = 10099, 11.4, 13.38, 0.015


def setup_path(which="latest"):
    sys.path.insert(0, LATEST if which == "latest" else GRIDCODE)


def load_ctx():
    from isochronetoolbox.scmd import Setup, context
    S = Setup(**{k: v for k, v in json.load(open(G + "/setup.json")).items()
                 if k in Setup.__dataclass_fields__})
    return S, context(S)


def window_points(ctx, S, cat):
    from isochronetoolbox.scmd.pipeline import _window
    p = ctx.points(cat)
    return p[_window(p, S.ks_window)]


def prep_target(pts, cdf=None):
    """(sorted points, their own ECDF = rank2d / n, sorted the same way)."""
    from isochronetoolbox.scmd.stats import rank2d
    if cdf is None:
        cdf = rank2d(pts) / len(pts)
    o = np.lexsort((pts[:, 0], pts[:, 1]))
    return pts[o], cdf[o]


def ks_many(sim, targets, q=95.0):
    """ks_2d(target, sim) for every target: the target's exact ECDF at its own points
    minus sim's ECDF interpolated (linear, 0 outside the hull) at the target's points.
    Returns a list of dicts with the q-th percentile, the 99th, the maximum and the mean."""
    from scipy.interpolate import LinearNDInterpolator
    from scipy.spatial import Delaunay
    from isochronetoolbox.scmd.stats import rank2d
    sim_cdf = rank2d(sim) / len(sim)
    f = LinearNDInterpolator(Delaunay(sim), sim_cdf, fill_value=0)
    res = []
    for p, c in targets:
        d = np.abs(c - f(p))
        res.append({"ks": float(np.percentile(d, q)), "p99": float(np.percentile(d, 99)),
                    "max": float(d.max()), "mean": float(d.mean())})
    return res


class Log:
    def __init__(self, name):
        self.fh = open(os.path.join(OUT, name + ".log"), "a")
        self.t0 = time.time()

    def __call__(self, *a):
        msg = f"[{time.time() - self.t0:7.1f} s] " + " ".join(str(x) for x in a)
        print(msg, flush=True)
        self.fh.write(msg + "\n")
        self.fh.flush()
