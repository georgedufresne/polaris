"""Helpers for the expected-ECDF audit (no repository file is modified).

cross_count / BigECDF: exact 2D dominance counts -- for each query point, the number
of reference ("big") points with colour strictly smaller AND magnitude strictly
smaller -- via a merge-sort tree over the colour-sorted reference set.  This is the
same quantity as rank2d(np.vstack([big, q]), mask)[~mask] (mask True for big) when
there are no ties between big and query coordinates, but O((Nb + Nq) log Nb) in
vectorised numpy instead of a Python recursion over a million points.
"""
import json
import sys

import numpy as np

sys.path.insert(0, '/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox')

G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
OUT = '/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/expected-ecdf'
RUN, AGE, MU, EVI = 10099, 11.4, 13.38, 0.015


class BigECDF:
    """Exact ECDF (strictly bluer AND brighter) of a fixed point set, queried anywhere."""

    def __init__(self, big: np.ndarray):
        big = np.asarray(big, float)
        order = np.argsort(big[:, 0], kind='stable')
        self.bc = big[order, 0]
        bm = big[order, 1]
        self.ub, rb = np.unique(bm, return_inverse=True)
        rb = rb.astype(np.int64).ravel()
        self.R = np.int64(len(self.ub) + 1)
        self.n = len(bm)
        idx = np.arange(self.n, dtype=np.int64)
        self.levels = []
        L = int(np.ceil(np.log2(max(self.n, 2)))) + 1
        for l in range(L):
            self.levels.append(np.sort((idx >> l) * self.R + rb))

    def count(self, q: np.ndarray) -> np.ndarray:
        q = np.asarray(q, float)
        p = np.searchsorted(self.bc, q[:, 0], side='left').astype(np.int64)   # colour strictly smaller
        rq = np.searchsorted(self.ub, q[:, 1], side='left').astype(np.int64)  # # distinct big mags < q_m
        cnt = np.zeros(len(q), dtype=np.int64)
        for l, keys in enumerate(self.levels):
            sel = ((p >> l) & 1).astype(bool)
            if not sel.any():
                continue
            start = (p[sel] >> (l + 1)) << (l + 1)
            b = start >> l
            cnt[sel] += np.searchsorted(keys, b * self.R + rq[sel], side='left') - start
        return cnt

    def __call__(self, q: np.ndarray) -> np.ndarray:
        return self.count(q) / self.n


def setup_context():
    from isochronetoolbox.scmd import Setup, context
    S = Setup(**{k: v for k, v in json.load(open(G + '/setup.json')).items()
                 if k in Setup.__dataclass_fields__})
    ctx = context(S)
    return S, ctx


def make_points(ctx, S, make, seed):
    from isochronetoolbox.scmd.pipeline import _window
    pts = ctx.points(make(np.random.RandomState(int(seed)), MU, EVI))
    return pts[_window(pts, S.ks_window)]
