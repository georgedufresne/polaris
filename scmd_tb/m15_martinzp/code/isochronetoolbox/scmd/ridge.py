"""Ridge line through a cluster CMD, and the cut around it.

G. Dufresne's NGC 6791 method, ported:

1. **warp** (optional) -- make the red stars temporarily brighter with a smooth
   monotonic function of colour (zero at ``c0``, PCHIP through the knots,
   capped), so that the sequence climbs steeply enough in magnitude for one
   colour per magnitude slice.  NGC 6791's subgiant branch even gets fainter to
   the red; 47 Tuc's brightens, but so little (0.1 in colour over 0.3 mag) that
   without a warp the DP cuts its corner, so the 47 Tuc preset uses one too;
2. **dynamic programme** -- in magnitude slices from faint to bright, the path
   minimising -log(KDE density) plus penalties on slope and curvature, with a
   per-slice limit on the colour step and an optional redward bias brighter
   than a magnitude (seam carving, Avidan & Shamir 2007);
3. **snap** the path to the local ridge: local PCA, weighted median along the
   normal, twice;
4. **unwarp** the ridge points;
5. **smooth** with a parametric spline and resample uniformly in arc length
   (as ``arc_length_continuous_ridge.csv``).

The sample is then the stars within a fixed distance of the ridge, measured
perpendicular to it with both axes normalised by their range (NGC 6791: 0.05).

For synthetic CMDs the same steps run on every sCMD (``RidgeCut``), so the
observed and the synthetic stars are selected in the same way.  That is ~2
million ridges per grid, hence the fast forms: the DP is vectorised over colour
(the same numbers as the loop, ``tests/scmd/test_ridge.py``), the snap over
ridge points, and the KDE is binned (``kde="exact"`` gives scikit-learn's).
The KDE bandwidth -- median distance to the 200th neighbour x 0.4 -- is fixed
from the observed stars and used for every sCMD.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence, Tuple

import numpy as np
from scipy.interpolate import PchipInterpolator, RegularGridInterpolator, splev, splprep
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree


def build_smooth_warp(c0: float, knots: Sequence[Tuple[float, float]], cap: float):
    """dG(colour): 0 at and blueward of ``c0``, PCHIP through ``knots``, never below ``cap``."""
    kx = [c0] + [t[0] for t in knots]
    ky = [0.0] + [max(cap, t[1]) for t in knots]
    warp = PchipInterpolator(kx, ky, extrapolate=True)

    def dG(x):
        return np.maximum(cap, np.minimum(0.0, warp(x)))
    return dG


def apply_warp(color, mag, dG):
    return mag + dG(color)


def inverse_warp_on_ridge(ridge_c, ridge_mw, dG):
    return ridge_mw - dG(ridge_c)


def bandwidth(colors, mags, bw_scale=0.4) -> float:
    """The DP's KDE bandwidth: median distance to the k-th neighbour (k = N/10 in 50..200) x ``bw_scale``."""
    X = np.column_stack([colors, mags])
    k = max(50, min(200, len(X) // 10))
    d, _ = cKDTree(X).query(X, k=k)
    return float(np.median(d[:, -1]) * bw_scale)


def turnoff_color(colors, mags, step=0.1, min_stars=50) -> float:
    """The bluest median colour of the 0.1-mag bins holding at least 50 stars."""
    edges = np.arange(mags.min(), mags.max() + step, step)
    k = np.digitize(mags, edges)
    return float(min(np.median(colors[k == b]) for b in np.unique(k) if np.sum(k == b) >= min_stars))


def _log_density(colors, mags, cols, mags_grid, bw, kde):
    if kde == "exact":
        from sklearn.neighbors import KernelDensity
        k = KernelDensity(bandwidth=bw, kernel="gaussian").fit(np.column_stack([colors, mags]))
        cc, gg = np.meshgrid(cols, mags_grid)
        return k.score_samples(np.column_stack([cc.ravel(), gg.ravel()])).reshape(len(mags_grid), len(cols))
    # binned: histogram at bw/3 per cell, Gaussian filter of 3 cells, linear interpolation to the grid
    dx = bw / 3.0
    xe = np.arange(cols.min() - 5 * bw, cols.max() + 5 * bw + dx, dx)
    ye = np.arange(mags_grid.min() - 5 * bw, mags_grid.max() + 5 * bw + dx, dx)
    h, _, _ = np.histogram2d(colors, mags, [xe, ye])
    dens = gaussian_filter(h, 3.0, mode="constant", truncate=5.0) / (len(colors) * dx * dx)
    f = RegularGridInterpolator((0.5 * (xe[1:] + xe[:-1]), 0.5 * (ye[1:] + ye[:-1])), dens)
    cc, gg = np.meshgrid(cols, mags_grid)
    d = f(np.column_stack([cc.ravel(), gg.ravel()])).reshape(len(mags_grid), len(cols))
    return np.log(np.maximum(d, 1e-12 * d.max()))


def ridge_via_dp(colors, mags, n_mag=220, n_col=260, col_pad=0.05, mag_pad=0.05, bw_scale=0.4,
                 dc_max=0.035, lam_slope=3.0, lam_curv=20.0, red_bias_after=None, lam_red=2.0,
                 bw=None, kde="exact"):
    """The dynamic-programming ridge: an (n, 2) array of (colour, mag), faint -> bright."""
    colors, mags = np.asarray(colors, float), np.asarray(mags, float)
    cmin, cmax = colors.min() - col_pad, colors.max() + col_pad
    gmax, gmin = mags.max() + mag_pad, mags.min() - mag_pad
    cols = np.linspace(cmin, cmax, n_col)
    mags_grid = np.linspace(gmax, gmin, n_mag)                 # slice i -> i+1 is brighter
    if bw is None:
        bw = bandwidth(colors, mags, bw_scale)
    cost = -_log_density(colors, mags, cols, mags_grid, bw, kde)

    bias = np.zeros_like(cost)
    if red_bias_after is not None:
        bright = mags_grid <= red_bias_after
        if bright.any():
            j0 = np.where(bright)[0][0]
            bias[j0:, :] -= (lam_red * np.maximum(0, cols - cols.mean()))[np.newaxis, :]

    J = np.full((n_mag, n_col), np.inf)
    prev = np.full((n_mag, n_col), -1, dtype=int)
    J[0, :] = cost[0, :]
    max_jump = max(1, int(np.ceil(dc_max / (cols[1] - cols[0]))))
    jj = np.arange(n_col)
    src = jj[np.newaxis, :] + np.arange(-max_jump, max_jump + 1)[:, np.newaxis]   # ascending source
    valid = (src >= 0) & (src < n_col)
    src = np.clip(src, 0, n_col - 1)
    dcol = cols[np.newaxis, :] - cols[src]
    for i in range(1, n_mag):
        slope = dcol / (mags_grid[i - 1] - mags_grid[i] + 1e-10)
        cand = J[i - 1][src] + cost[i][np.newaxis, :] + bias[i][np.newaxis, :] + lam_slope * slope ** 2
        if i >= 2:
            curv = cols[np.newaxis, :] - 2 * cols[src] + cols[prev[i - 1][src]]
            cand = cand + lam_curv * curv ** 2
        cand[~valid] = np.inf
        kb = np.argmin(cand, axis=0)
        J[i] = cand[kb, jj]
        prev[i] = src[kb, jj]
    j = int(np.argmin(J[-1, :]))
    path = [j]
    for i in range(n_mag - 1, 0, -1):
        j = prev[i, j]
        if j < 0 or j >= n_col:
            break
        path.append(j)
    path = path[::-1]
    return np.column_stack([cols[np.array(path)], mags_grid[:len(path)]])


def snap_to_local_ridge(ridge, colors, mags, k_neighbors=80, n_iter=2, tree=None):
    """Move each ridge point to the weighted median of its neighbours along the local normal."""
    X = np.column_stack([colors, mags])
    tree = cKDTree(X) if tree is None else tree
    snapped = np.array(ridge, float)
    for _ in range(n_iter):
        dist, idx = tree.query(snapped, k=k_neighbors)
        nb = X[idx]                                                    # (n, k, 2)
        dev = nb - nb.mean(axis=1, keepdims=True)
        _, vec = np.linalg.eigh(np.einsum("nki,nkj->nij", dev, dev))
        normal = vec[:, :, 0]                                          # minor axis, as PCA's 2nd component
        proj = np.einsum("nki,ni->nk", nb - snapped[:, np.newaxis, :], normal)
        w = np.exp(-dist / dist.mean(axis=1, keepdims=True))
        o = np.argsort(proj, axis=1)
        cs = np.cumsum(np.take_along_axis(w, o, axis=1), axis=1)
        m = np.minimum((cs < 0.5 * cs[:, -1:]).sum(axis=1), k_neighbors - 1)
        off = np.take_along_axis(proj, np.take_along_axis(o, m[:, np.newaxis], axis=1), axis=1)
        snapped = snapped + off * normal
    return snapped


def _arc_resample(c, m, n):
    s = np.r_[0, np.cumsum(np.hypot(np.diff(c), np.diff(m)))]
    u = np.linspace(0, s[-1], n)
    return np.column_stack([np.interp(u, s, c), np.interp(u, s, m)])


def smooth_ridge(ridge, smoothing=3e-5, n=1000):
    """Parametric smoothing spline (s = ``smoothing`` x points) through the ordered
    ridge, resampled at ``n`` points uniform in arc length."""
    keep = np.r_[True, np.any(np.diff(ridge, axis=0) != 0, axis=1)]   # splprep refuses repeats
    r = ridge[keep]
    tck, _ = splprep([r[:, 0], r[:, 1]], s=smoothing * len(r), k=3)
    c, m = splev(np.linspace(0, 1, 20 * n), tck)
    return _arc_resample(c, m, n)


@dataclass
class Tracer:
    """The settings of one ridge: warp, DP, snap, smoothing.  Defaults are the NGC 6791 run's."""
    warp: Optional[Tuple[float, Tuple[Tuple[float, float], ...], float]] = None   # (c0, knots, cap)
    warp_from_turnoff: bool = False      # warp colours are offsets from the CMD's own turnoff colour
    n_mag: int = 120
    n_col: int = 240
    dc_max: float = 0.035
    lam_slope: float = 3.0
    lam_curv: float = 15.0
    red_bias_after: Optional[float] = None
    lam_red: float = 1.5
    snap_k: int = 100
    snap_iter: int = 2
    smoothing: float = 3e-5
    n: int = 1000
    kde: str = "binned"

    def __call__(self, colors, mags, bw=None) -> np.ndarray:
        """Warp -> DP -> snap -> unwarp -> smooth: the ridge as (n, 2) (colour, mag), faint -> bright."""
        colors, mags = np.asarray(colors, float), np.asarray(mags, float)
        if self.warp:
            c0, knots, cap = self.warp
            if self.warp_from_turnoff:
                t = turnoff_color(colors, mags)
                c0, knots = c0 + t, tuple((t + kc, kg) for kc, kg in knots)
            dG = build_smooth_warp(c0, knots, cap)
        else:
            dG = lambda x: np.zeros_like(np.asarray(x, float))        # noqa: E731
        mw = apply_warp(colors, mags, dG)
        dp = ridge_via_dp(colors, mw, n_mag=self.n_mag, n_col=self.n_col, dc_max=self.dc_max,
                          lam_slope=self.lam_slope, lam_curv=self.lam_curv,
                          red_bias_after=self.red_bias_after, lam_red=self.lam_red, bw=bw, kde=self.kde)
        sn = snap_to_local_ridge(dp, colors, mw, k_neighbors=self.snap_k, n_iter=self.snap_iter)
        ridge = np.column_stack([sn[:, 0], inverse_warp_on_ridge(sn[:, 0], sn[:, 1], dG)])
        return smooth_ridge(ridge, self.smoothing, self.n)


class Ridge:
    """A ridge line and the perpendicular distance to it in normalised CMD units."""

    def __init__(self, color, mag, scales: Tuple[float, float], n: int = 4000):
        self.color, self.mag = np.asarray(color, float), np.asarray(mag, float)
        x, y = self.color / scales[0], self.mag / scales[1]
        s = np.r_[0, np.cumsum(np.hypot(np.diff(x), np.diff(y)))]
        u = np.linspace(0, s[-1], n)
        self.tree = cKDTree(np.column_stack([np.interp(u, s, x), np.interp(u, s, y)]))
        self.scales = scales

    def distance(self, color, mag) -> np.ndarray:
        d, _ = self.tree.query(np.column_stack([np.asarray(color) / self.scales[0],
                                                np.asarray(mag) / self.scales[1]]))
        return d


@dataclass
class RidgeCut:
    """The fit-star selection, for observed and synthetic CMDs alike: a ridge traced
    through the ``n_trace`` stars in the colour box and the magnitude window widened
    by ``margin`` (the DP's first and last slices are unreliable, so the ridge is
    traced past the window and trimmed to it), then the stars in the window and box
    within ``width`` of the ridge (axes divided by ``scales``).  Stars inside any
    ``exclude`` box ((colour range), (mag range)) are left out of both, as the NGC
    6791 run left out the blue stragglers' corner before tracing."""
    tracer: Tracer
    mag_range: Tuple[float, float]
    color_range: Tuple[float, float]
    width: float
    margin: float = 0.3
    exclude: Tuple = ()
    scales: Tuple[float, float] = (1.0, 1.0)
    bw: Optional[float] = None
    n_trace: int = 0
    ridge: Optional[Ridge] = field(default=None, repr=False)

    def in_box(self, color, mag, margin: float = 0.0) -> np.ndarray:
        ok = ((mag > self.mag_range[0] - margin) & (mag < self.mag_range[1] + margin)
              & (color > self.color_range[0]) & (color < self.color_range[1]))
        for (c0, c1), (m0, m1) in self.exclude:
            ok &= ~((color > c0) & (color < c1) & (mag > m0) & (mag < m1))
        return ok

    def in_trace_box(self, color, mag) -> np.ndarray:
        return self.in_box(color, mag, self.margin)

    def trace(self, color, mag) -> Ridge:
        r = self.tracer(color, mag, self.bw)
        r = r[(r[:, 1] >= self.mag_range[0]) & (r[:, 1] <= self.mag_range[1])]
        return Ridge(r[:, 0], r[:, 1], self.scales)

    def select(self, color, mag, ridge: Ridge) -> np.ndarray:
        keep = self.in_box(color, mag)
        keep[keep] = ridge.distance(color[keep], mag[keep]) <= self.width
        return keep

    @classmethod
    def from_observed(cls, tracer: Tracer, color, mag, mag_range, color_range, width, margin=0.3,
                      scales: Optional[Tuple[float, float]] = None, exclude=()) -> Tuple["RidgeCut", np.ndarray]:
        """The cut fixed from the observed stars -- normalisation (their range in the window
        and box unless ``scales`` is given), bandwidth, star count -- and the mask of the kept stars."""
        color, mag = np.asarray(color, float), np.asarray(mag, float)
        rc = cls(tracer, tuple(mag_range), tuple(color_range), width, margin,
                 tuple((tuple(c), tuple(m)) for c, m in exclude))
        box = rc.in_box(color, mag)
        rc.scales = tuple(scales) if scales else (float(np.ptp(color[box])), float(np.ptp(mag[box])))
        tb = rc.in_trace_box(color, mag)
        rc.bw = bandwidth(color[tb], mag[tb])
        rc.n_trace = int(tb.sum())
        rc.ridge = rc.trace(color[tb], mag[tb])
        return rc, rc.select(color, mag, rc.ridge)
