"""Comparing a synthetic CMD with the observed one: KDE Jensen-Shannon and 2D KS.

Condensed from the stage-2 (``kde_code_for_grid_best_match_chaboyer_test.py``)
and stage-3 (``2d_ks_comparison_grid_chaboyer_kde.py``) scripts.

:class:`KDE` reproduces ``KDEpy.FFTKDE(kernel='gaussian', bw, norm=1)`` evaluated
on a fixed 1000 x 1000 grid, to the bit: same linear binning (KDEpy's own
Cython routine), same kernel weights (KDEpy's own kernel), and the same FFT
convolution as ``scipy.signal.convolve(mode='same')``.  It is ~6x faster
because the grid bookkeeping (``np.unique`` over a million points, twice per
call), the 0.7-million-point kernel and its FFT are computed once instead of
on every call, and the FFTs skip the rows that are only padding.
"""

from __future__ import annotations

import numpy as np
from KDEpy import _cutils
from KDEpy.kernel_funcs import gaussian
from KDEpy.utils import cartesian
from scipy.fft import fft, ifft, irfft, next_fast_len, rfft, rfftn
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial.distance import jensenshannon

AXES = [0, 1]


class KDE:
    """Gaussian KDE (L1 norm, as in the original) on a fixed equidistant grid."""

    def __init__(self, color_range, mag_range, n: int = 1000, bw: float = 0.15,
                 eps: float = 1e-6):
        cg = np.linspace(color_range[0] - eps, color_range[1] + eps, n)
        mg = np.linspace(mag_range[0] - eps, mag_range[1] + eps, n)
        cm, mm = np.meshgrid(cg, mg)
        pts = np.column_stack([cm.ravel(), mm.ravel()])
        pts = pts[np.lexsort((pts[:, 1], pts[:, 0]))]
        self.lo, self.hi = np.min(pts, axis=0), np.max(pts, axis=0)
        self.num = np.array([len(np.unique(pts[:, i])) for i in range(2)], dtype="long")
        self.dx = (self.hi - self.lo) / (self.num - 1)
        # kernel exactly as FFTKDE.evaluate builds it
        L = np.minimum(np.floor(gaussian.practical_support(bw) / self.dx), self.num)
        grids = [np.linspace(-d * l, d * l, int(l * 2 + 1)) for d, l in zip(self.dx, L)]
        w = gaussian(cartesian(grids), bw=bw, norm=1).reshape(*[int(l * 2 + 1) for l in L])
        # scipy.signal.fftconvolve(binned, w, mode='same') with w's FFT cached
        full = [int(a + b - 1) for a, b in zip(self.num, w.shape)]
        self.fshape = [next_fast_len(s, True) for s in full]
        self.W = rfftn(w, self.fshape, axes=AXES)
        start = (np.array(full) - self.num) // 2
        self.same = tuple(slice(int(s), int(s + k)) for s, k in zip(start, self.num))
        # irfftn's 1/N, as pocketfft computes it (long double, applied once at the end)
        self.fct = np.float64(np.longdouble(1) / np.longdouble(self.fshape[0] * self.fshape[1]))

    def pdf(self, data: np.ndarray) -> np.ndarray:
        """Density on the grid (flattened, color-major), as FFTKDE returns it."""
        data = np.asarray_chkfinite(data, dtype=float)
        if not ((self.lo < data.min(axis=0)).all() and (self.hi > data.max(axis=0)).all()):
            raise ValueError("Every data point must be inside of the grid.")
        binned = np.zeros(int(self.num[0] * self.num[1]), dtype=float)
        binned = _cutils.iterate_data_2D((data - self.lo) / self.dx, binned, self.num, len(binned))
        binned = np.asarray_chkfinite(binned, dtype=float) / data.shape[0]
        # rfftn/irfftn one axis at a time, skipping the rows that are zero on the
        # way in (padding) or discarded on the way out (outside 'same'): the
        # remaining arithmetic is pocketfft's, so the result is unchanged
        F0, F1 = self.fshape
        sp = fft(rfft(binned.reshape(*tuple(self.num)), n=F1, axis=1), n=F0, axis=0)
        rows = ifft(sp * self.W, axis=0, norm="forward")[self.same[0]]
        conv = irfft(rows, n=F1, axis=1, norm="forward")[:, self.same[1]]
        return (conv * self.fct).ravel() + np.finfo(float).eps


def normalise(pdf: np.ndarray) -> np.ndarray:
    pdf = np.maximum(pdf, 1e-10)
    return pdf / np.sum(pdf)


def js_masked(p_obs: np.ndarray, p_sim: np.ndarray) -> float:
    """JS distance over the cells where both normalised densities exceed 1e-8."""
    mask = (p_obs > 1e-8) & (p_sim > 1e-8)
    if np.sum(mask) < 100:
        return np.inf
    js = jensenshannon(p_obs[mask], p_sim[mask])
    return np.inf if (np.isnan(js) or np.isinf(js)) else js


def js_fixed(kde: KDE, p_obs: np.ndarray, sim: np.ndarray) -> float:
    """Stage 2: JS distance on the fixed grid against a precomputed observed pdf."""
    sim = sim[np.isfinite(sim).all(axis=1)]
    if len(sim) < 10:
        return np.inf
    try:
        return js_masked(p_obs, normalise(kde.pdf(sim)))
    except Exception:
        return np.inf


def js_adaptive(obs: np.ndarray, sim: np.ndarray, n: int = 1000, bw: float = 0.15,
                eps: float = 1e-6) -> float:
    """Stage 3: JS distance on a 1000x1000 grid spanning both data sets."""
    both = np.vstack([obs, sim])
    obs = obs[np.isfinite(obs).all(axis=1)]
    sim = sim[np.isfinite(sim).all(axis=1)]
    if len(obs) < 10 or len(sim) < 10:
        return np.inf
    try:
        kde = KDE((both[:, 0].min(), both[:, 0].max()), (both[:, 1].min(), both[:, 1].max()), n, bw, eps)
        return js_masked(normalise(kde.pdf(obs)), normalise(kde.pdf(sim)))
    except Exception:
        return np.inf


def rank2d(points: np.ndarray, mask=None) -> np.ndarray:
    """Number of points strictly below-left of each point (divide and conquer)."""
    N = points.shape[0]
    N2 = N // 2
    if N == 1:
        return np.array([0])
    idx = np.argpartition(points[:, 0], N2)
    idxA = np.zeros(N, dtype=bool)
    idxA[idx[:N2]] = True
    if mask is not None:
        NAm = np.sum(idxA & mask)
        reduced = np.vstack((points[idxA & mask], points[~idxA & ~mask]))
    else:
        NAm = np.sum(idxA)
        reduced = np.vstack((points[idxA], points[~idxA]))
    count = np.zeros(reduced.shape[0], dtype=bool)
    count[:NAm] = True
    idxY = np.argsort(reduced[:, 1])
    idxYr = np.zeros_like(idxY)
    idxYr[idxY] = np.arange(idxY.shape[0])
    numA = np.cumsum(count[idxY])[idxYr]
    rank = np.zeros(N, dtype=int)
    if mask is not None:
        rank[idxA] = rank2d(points[idxA], mask[idxA])
        rank[~idxA] = rank2d(points[~idxA], mask[~idxA])
        rank[~idxA & ~mask] += numA[NAm:]
    else:
        rank[idxA] = rank2d(points[idxA])
        rank[~idxA] = rank2d(points[~idxA])
        rank[~idxA] += numA[NAm:]
    return rank


def ks_2d(obs: np.ndarray, sim: np.ndarray, obs_cdf: np.ndarray = None, q: float = 95.0):
    """``q``-th percentile (default the 95th, as the NGC 6791 pipeline; M. Ying used
    the 99th) of |ECDF_obs - ECDF_sim| at the observed points.

    The simulated ECDF is interpolated linearly (0 outside its hull).  Returns
    (statistic, differences).  ``obs_cdf`` may be passed to reuse the observed
    ranks across many sCMDs.
    """
    if obs_cdf is None:
        obs_cdf = rank2d(obs) / len(obs)
    sim_cdf = rank2d(sim) / len(sim)
    d = np.abs(obs_cdf - LinearNDInterpolator(sim, sim_cdf, fill_value=0)(obs))
    return np.percentile(d, q), d
