"""Synthetic CMDs with artificial-star photometry (HST clusters, e.g. 47 Tuc).

The NGC 6791 model (``model.py``) injects empirical scatter perpendicular to
the isochrone.  For HST fields the photometric errors and the completeness
come from artificial-star tests instead, and both depend on how crowded the
star's position is, i.e. on its distance from the cluster centre.  This module
follows the error model of M. Ying's sCMD code (``PhotoPropI.f90``,
``RadialDensity.f90``, ``MakeTestCMD.f90``):

* the artificial stars are binned, separately in F606W and F814W, in 10
  radial x 12 magnitude bins (6 for the short and 6 for the long exposures);
  each bin keeps its completeness and its sorted list of (input - output)
  magnitudes for the stars that were recovered;
* each synthetic star gets a distance from the centre drawn from the observed
  fit stars' radial distribution, and in each band a completeness test and an
  error drawn from the matching bin.

What differs from that code, and why:

* magnitudes are placed in the artificial-star frame with the photometric zero
  points (``mag_inst = mag - zp``) rather than by matching the subgiant branch
  of each isochrone, because here the grid applies the distance modulus and
  reddening to the sCMD (as for NGC 6791), so its stars have real apparent
  magnitudes;
* the error is applied as ``observed = true - (input - output)``; the Fortran
  adds (input - output), which flips any asymmetry of the error distribution
  (blends come out fainter instead of brighter).  ``martin_sign=True``
  reproduces it;
* the sCMD gets the observed fit stars' selection -- the F606W window and a
  0.08-mag cut in both bands, around the sCMD's own isochrone as in
  MakeTestCMD (the observed stars were cut around their fiducial line) -- and
  as many stars as the observed sample has.  Cutting the sCMD around the
  *observed* fiducial instead (``cut="fiducial"``) is a trap with a fixed star
  count: a misplaced isochrone is filled up with its error-tail stars.

The statistics afterwards are the NGC 6791 ones (``stats.py``).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline, interp1d

F32 = np.float32


def centre(x: np.ndarray, y: np.ndarray, guess=(3000.0, 3000.0), radius=200.0) -> Tuple[float, float]:
    """Mean position of the stars within ``radius`` px of ``guess`` (as in the Fortran)."""
    x, y = x.astype(F32), y.astype(F32)
    near = (x - F32(guess[0])) ** 2 + (y - F32(guess[1])) ** 2 <= F32(radius) ** 2
    n = F32(near.sum())
    # summed one star at a time in single precision, as the Fortran loop does
    return float(np.cumsum(x[near], dtype=F32)[-1] / n), float(np.cumsum(y[near], dtype=F32)[-1] / n)


@dataclass
class BandErrors:
    """Completeness and error distributions of one band on the (radius, mag) grid."""
    rad_edges: np.ndarray            # outer edge of each radial bin (Nrad)
    mag_edges: np.ndarray            # faint edge of each magnitude bin (Nmag), instrumental
    completeness: np.ndarray         # (Nrad, Nmag)
    errors: Dict[Tuple[int, int], np.ndarray]   # sorted (input - output) per bin

    def bins(self, radius: np.ndarray, mag_inst: np.ndarray):
        """Bin indices as MakeTestCMD's findIndxForErrors assigns them (0-based)."""
        nr, nm = len(self.rad_edges), len(self.mag_edges)
        ir = np.searchsorted(self.rad_edges[: nr - 1], radius, side="left")      # first edge >= r
        im = np.searchsorted(self.mag_edges[: nm - 1], mag_inst, side="left")
        return np.minimum(ir, nr - 1), np.minimum(im, nm - 1)

    def draw(self, rs: np.random.RandomState, radius: np.ndarray, mag_inst: np.ndarray):
        """(recovered, error) for each star: a completeness test, then an error
        from the bin's distribution by inverse CDF (cdf_k = k/n, k = 1..n)."""
        ir, im = self.bins(radius, mag_inst)
        ok = rs.random_sample(len(radius)) < self.completeness[ir, im]
        err = np.zeros(len(radius))
        for key in set(zip(ir[ok].tolist(), im[ok].tolist())):
            sel = np.where(ok & (ir == key[0]) & (im == key[1]))[0]
            e = self.errors[key]
            cdf = np.arange(1, len(e) + 1) / len(e)
            u = rs.random_sample(len(sel))
            bad = u < cdf[0]
            while bad.any():                      # the Fortran redraws below the first step
                u[bad] = rs.random_sample(bad.sum())
                bad = u < cdf[0]
            err[sel] = np.interp(u, cdf, e)
        return ok, err


def bin_artificial_stars(x_in, y_in, m_in, x_out, y_out, m_out, mag_min, mag_max, short_long,
                         n_rad=10, n_mag=12, guess=(3000.0, 3000.0), guess_radius=200.0,
                         pix_miss=0.5, mag_miss=0.75) -> BandErrors:
    """PhotoPropI.f90 for one band, in single precision as the Fortran computes it.

    A star is recovered if it is found within 0.5 px and 0.75 mag of its input
    (Anderson et al. 2008).  The radial bins start 1.5 widths out and span the
    input field; the magnitude bins split at the short/long exposure boundary.
    """
    sel = (m_in <= F32(mag_max)) & (m_in >= F32(mag_min))
    x_in, y_in, m_in = x_in[sel], y_in[sel], m_in[sel]
    x_out, y_out, m_out = x_out[sel], y_out[sel], m_out[sel]
    rec = (np.sqrt((x_in - x_out) ** 2 + (y_in - y_out) ** 2) < F32(pix_miss)) & \
          (np.abs(m_in - m_out) < F32(mag_miss))
    xc, yc = centre(x_in, y_in, guess, guess_radius)
    dist = np.sqrt((x_in - F32(xc)) ** 2 + (y_in - F32(yc)) ** 2).astype(F32)

    maxd = dist.max()
    w = maxd / F32(n_rad)
    rad = np.zeros(n_rad, dtype=F32)
    rad[0], rad[-1] = w * F32(1.5), maxd
    w = (maxd - F32(3.0) * w) / F32(n_rad)
    for i in range(1, n_rad - 1):
        rad[i] = rad[i - 1] + w
    vmin, vmax, vcut = m_in.min(), m_in.max(), F32(short_long)
    w1 = (vcut - vmin) / F32(n_mag) * F32(2.0)
    w2 = (vmax - vcut) / F32(n_mag) * F32(2.0)
    half = n_mag // 2
    mag = np.array([vmin + F32(i) * w1 for i in range(1, half + 1)] +
                   [vcut + F32(i - half) * w2 for i in range(half + 1, n_mag + 1)], dtype=F32)

    ir = np.where(dist <= rad[0], 0,
                  np.where(dist > rad[n_rad - 2], n_rad - 1,
                           ((dist - rad[0]) / w).astype(np.int32) + 1))
    im = np.where(m_in <= vcut, ((m_in - vmin) / w1).astype(np.int32),
                  np.minimum(((m_in - vcut) / w2).astype(np.int32) + half, n_mag - 1))
    comp = np.zeros((n_rad, n_mag))
    errors = {}
    for i in range(n_rad):
        for j in range(n_mag):
            b = (ir == i) & (im == j)
            good = b & rec
            ng, nb = int(good.sum()), int((b & ~rec).sum())
            comp[i, j] = F32(ng) / F32(ng + nb) if ng + nb else 0.0
            errors[(i, j)] = np.sort((m_in[good] - m_out[good]).astype(F32)).astype(float)
    return BandErrors(rad.astype(float), mag.astype(float), comp, errors)


class ArtificialStars:
    """Both bands' error model plus the radial distribution of the observed stars."""

    def __init__(self, ast_csv: str, fit_stars: pd.DataFrame, zp: Dict[str, float],
                 window: Dict[str, Tuple[float, float]], short_long: Dict[str, float],
                 n_rad=10, n_mag=12, guess=(3000.0, 3000.0), guess_radius=200.0):
        a = pd.read_csv(ast_csv, dtype=F32)
        self.zp = zp
        self.band = {}
        for b, (ci, co) in {"F606W": ("InputF606W", "OutputF606W"),
                            "F814W": ("InputF814W", "OutputF814W")}.items():
            lo, hi = window[b]
            self.band[b] = bin_artificial_stars(
                a["InputX"].values, a["InputY"].values, a[ci].values,
                a["OutputX"].values, a["OutputY"].values, a[co].values,
                lo - zp[b], hi - zp[b], short_long[b] - zp[b], n_rad, n_mag, guess, guess_radius)
        # RadialDensity.f90: distances of the observed fit stars from their own centre
        x, y = fit_stars["x"].values, fit_stars["y"].values
        xc, yc = centre(x, y, guess, guess_radius)
        self.dist = np.sort(np.sqrt((x.astype(F32) - F32(xc)) ** 2 + (y.astype(F32) - F32(yc)) ** 2)).astype(float)
        self.cdf = np.arange(1, len(self.dist) + 1) / len(self.dist)

    def radii(self, rs: np.random.RandomState, n: int) -> np.ndarray:
        """Distances drawn from the observed radial distribution (inverse CDF)."""
        u = rs.random_sample(n)
        bad = u < self.cdf[0]
        while bad.any():
            u[bad] = rs.random_sample(bad.sum())
            bad = u < self.cdf[0]
        return np.interp(u, self.cdf, self.dist)


class Fiducial:
    """The observed fiducial line and the 0.08-mag band around it that selects the fit stars."""

    def __init__(self, path: str, width: float = 0.08):
        f = pd.read_csv(path)
        vp, ip = f["v"].values, f["v"].values - f["vi"].values
        vp = vp[::-1] if vp[0] > vp[-1] else vp
        ip = ip[::-1] if ip[0] > ip[-1] else ip
        self.v_of_i = CubicSpline(ip, vp, extrapolate=True)
        self.i_of_v = CubicSpline(vp, ip, extrapolate=True)
        self.width = width

    def keep(self, v: np.ndarray, i: np.ndarray) -> np.ndarray:
        dv, di = self.v_of_i(i), self.i_of_v(v)
        w = self.width
        return (dv + w > v) & (v > dv - w) & (di + w > i) & (i > di - w)


class Ridge:
    """A ridge line in the CMD and the perpendicular distance to it.

    The NGC 6791 technique: distances are measured in normalised CMD units
    (colour / ``scales[0]``, magnitude / ``scales[1]``), perpendicular to the
    ridge, and the sample is the stars within a fixed normalised distance.
    """

    def __init__(self, color: np.ndarray, mag: np.ndarray, scales: Tuple[float, float],
                 smooth: Optional[float] = None, n: int = 4000):
        x, y = np.asarray(color, float) / scales[0], np.asarray(mag, float) / scales[1]
        if smooth is None:          # an isochrone: the polyline itself, densified
            s = np.r_[0, np.cumsum(np.hypot(np.diff(x), np.diff(y)))]
            u = np.linspace(0, s[-1], n)
            self.x, self.y = np.interp(u, s, x), np.interp(u, s, y)
        else:                       # data: a smoothing parametric spline through the points
            from scipy.interpolate import splev, splprep
            tck, _ = splprep([x, y], s=smooth, k=3)
            self.x, self.y = splev(np.linspace(0, 1, n), tck)
        from scipy.spatial import cKDTree
        self.tree = cKDTree(np.column_stack([self.x, self.y]))
        self.scales = scales
        self.points = (np.asarray(color, float), np.asarray(mag, float))

    def distance(self, color: np.ndarray, mag: np.ndarray) -> np.ndarray:
        d, _ = self.tree.query(np.column_stack([np.asarray(color) / self.scales[0],
                                                np.asarray(mag) / self.scales[1]]))
        return d

    @classmethod
    def from_density(cls, color: np.ndarray, mag: np.ndarray, scales: Tuple[float, float],
                     mag_range: Tuple[float, float], color_range: Tuple[float, float] = (0.35, 0.95),
                     bins: Tuple[int, int] = (180, 400), sigma: float = 2.0) -> "Ridge":
        """A ridge line through an observed cluster, after the NGC 6791 method.

        That method follows the CMD's KDE density with a seam-carving dynamic
        programme (penalising changes of slope and curvature, with a temporary
        local function for the subgiant bend).  Stripped down: the cheapest path
        through the smoothed density, in normalised units, from the densest
        point of the faint end to the densest point of the bright end.  A path
        over the whole grid (not one step per magnitude row) turns along the
        subgiant branch by itself; a smoothing spline then removes the pixel steps.
        """
        from scipy.ndimage import gaussian_filter
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import dijkstra
        s0, s1 = scales
        xe = np.linspace(color_range[0] / s0, color_range[1] / s0, bins[0] + 1)
        ye = np.linspace(mag_range[0] / s1, mag_range[1] / s1, bins[1] + 1)
        h, _, _ = np.histogram2d(np.asarray(color) / s0, np.asarray(mag) / s1, bins=[xe, ye])
        dens = gaussian_filter(h, sigma)
        cost = 1.0 / (dens + 1e-3 * dens.max())
        nx, ny = dens.shape
        dx, dy = xe[1] - xe[0], ye[1] - ye[0]
        idx = np.arange(nx * ny).reshape(nx, ny)
        rows, cols, w = [], [], []
        for ox, oy in ((1, 0), (0, 1), (1, 1), (1, -1)):
            a = idx[max(0, -ox):nx - max(0, ox), max(0, -oy):ny - max(0, oy)]
            b = idx[max(0, ox):nx - max(0, -ox) or None, max(0, oy):ny - max(0, -oy) or None]
            step = np.hypot(ox * dx, oy * dy)
            c = 0.5 * (cost.ravel()[a.ravel()] + cost.ravel()[b.ravel()]) * step
            rows += [a.ravel(), b.ravel()]
            cols += [b.ravel(), a.ravel()]
            w += [c, c]
        g = coo_matrix((np.concatenate(w), (np.concatenate(rows), np.concatenate(cols))), shape=(nx * ny,) * 2).tocsr()
        start = idx[np.argmax(dens[:, -3:].sum(axis=1)), ny - 2]     # faint end
        end = idx[np.argmax(dens[:, :3].sum(axis=1)), 1]               # bright end
        _, pred = dijkstra(g, indices=start, return_predecessors=True)
        path, k = [], end
        while k != start and k >= 0:
            path.append(k)
            k = pred[k]
        path = np.array(path[::-1])
        px = 0.5 * (xe[:-1] + xe[1:])[path // ny] * s0
        py = 0.5 * (ye[:-1] + ye[1:])[path % ny] * s1
        return cls(px[::-1], py[::-1], scales, smooth=len(px) * 3e-5)

    @classmethod
    def from_data(cls, color: np.ndarray, mag: np.ndarray, scales: Tuple[float, float],
                  ms_top: float, rgb_bottom: float, mag_range: Tuple[float, float],
                  ms_bin: float = 0.05, rgb_bin: float = 0.08, sgb_bin: float = 0.01) -> "Ridge":
        """A fiducial line through an observed cluster: sigma-clipped median colour
        in magnitude bins on the main sequence (fainter than ``ms_top``) and the
        giant branch (brighter than ``rgb_bottom``), median magnitude in colour
        bins along the flat subgiant branch in between."""
        def med(x, y, edges):
            out = []
            for a, b in zip(edges[:-1], edges[1:]):
                yy = y[(x >= a) & (x < b)]
                for _ in range(3):
                    if len(yy) < 20:
                        break
                    m = np.median(yy)
                    yy = yy[np.abs(yy - m) < 3 * max(1.4826 * np.median(np.abs(yy - m)), 1e-3)]
                if len(yy) >= 20:
                    out.append((0.5 * (a + b), np.median(yy)))
            return np.array(out)
        lo, hi = mag_range
        ms = med(mag, color, np.arange(ms_top, hi + ms_bin / 2, ms_bin))[:, ::-1]
        rgb = med(mag, color, np.arange(lo, rgb_bottom + rgb_bin / 2, rgb_bin))[:, ::-1]
        ms, rgb = ms[np.argsort(-ms[:, 1])], rgb[np.argsort(-rgb[:, 1])]
        band = (mag > rgb_bottom - 0.05) & (mag < ms_top + 0.05)
        to = ms[ms[:, 1] < ms_top + 0.15, 0].min()
        sgb = med(color[band], mag[band], np.arange(to + 2 * sgb_bin, rgb[-1, 0] - sgb_bin / 2, sgb_bin))
        pts = np.vstack([ms, sgb, rgb])
        return cls(pts[:, 0], pts[:, 1], scales, smooth=len(pts) * 2e-5)


def lininterp(x: np.ndarray, y: np.ndarray, xv: np.ndarray) -> np.ndarray:
    """MakeTestCMD's ``lininterp`` + ``locate``, vectorised.

    Linear interpolation with a bisection search whose direction is set by the
    end points, clamped to the end values outside them.  On a non-monotonic x
    (V along the subgiant branch) it picks the crossing the bisection lands
    on, exactly as the Fortran does.
    """
    n = len(x)
    xv = np.asarray(xv, dtype=float)
    asc = x[-1] >= x[0]
    jl, ju = np.zeros(len(xv), dtype=int), np.full(len(xv), n + 1)
    while True:
        act = ju - jl > 1
        if not act.any():
            break
        jm = (ju + jl) // 2
        go = act & ((xv >= x[np.clip(jm - 1, 0, n - 1)]) == asc)
        jl = np.where(go, jm, jl)
        ju = np.where(act & ~go, jm, ju)
    i = jl.copy()                                    # 1-based, as in the Fortran
    i[xv == x[0]], i[xv == x[-1]] = 1, n - 1
    if asc:
        i = np.clip(i, 1, n - 1)
        val = (y[i] - y[i - 1]) / (x[i] - x[i - 1]) * (xv - x[i - 1]) + y[i - 1]
        return np.where(xv <= x[0], y[0], np.where(xv >= x[-1], y[-1], val))
    i = np.clip(i, 2, n)
    val = (y[i - 1] - y[i - 2]) / (x[i - 1] - x[i - 2]) * (xv - x[i - 2]) + y[i - 2]
    return np.where(xv <= x[-1], y[-1], np.where(xv >= x[0], y[0], val))


def make_scmd_ast(rs: np.random.RandomState, mass: np.ndarray, phot: Dict[str, np.ndarray],
                  ast: ArtificialStars, fid: Optional[Fiducial], n_target: int,
                  v_window: Tuple[float, float], pad: float = 1.0, imf_slope: float = -0.84,
                  binary_fraction: float = 0.0, martin_sign: bool = False, cut: str = "isochrone",
                  width: float = 0.08, max_factor: float = 50.0,
                  ridge_scales: Tuple[float, float] = (0.45, 4.0)) -> Optional[np.ndarray]:
    """One synthetic CMD with artificial-star errors: an (n, 2) array of observed
    (F606W, F814W), n = ``n_target`` unless the cell is hopeless; None if the
    isochrone misses the window.

    Masses come from the power-law mass function between the lowest and highest
    isochrone mass with F606W within ``pad`` mag of the window, so stars can
    scatter in from outside it.  Each star gets a distance from the centre and,
    in each band, a completeness test and an artificial-star error.  Then the
    selection of the observed fit stars: the F606W window, and a cut around
    the sCMD's own isochrone -- ``cut="ridge"``: within ``width`` in normalised
    perpendicular distance, as the observed stars are cut around their own
    ridge line (``Ridge.from_data``); ``cut="isochrone"``: 0.08 mag in both
    bands as in MakeTestCMD -- or around the observed fiducial (``"fiducial"``).  Stars are
    drawn in batches until ``n_target`` survive, up to ``max_factor`` x
    ``n_target`` draws.
    """
    V, I = phot["F606W"], phot["F814W"]
    o = np.argsort(mass)
    m, Vs, Is = mass[o], V[o], I[o]
    lo, hi = v_window
    win = np.where((Vs >= lo - pad) & (Vs <= hi + pad))[0]
    if not len(win):
        return None
    a, b = win[0], win[-1]                           # contiguous range, as indxmin:indxmax
    m_lo, m_hi = m[a:b + 1].min(), m[a:b + 1].max()
    # MakeTestCMD passes npts = indxmax - indxmin, one short of the slice
    ri, rv = Is[a:b], Vs[a:b]
    ridge = Ridge(rv - ri, rv, ridge_scales) if cut == "ridge" else None
    fV = interp1d(m, Vs, kind="linear", fill_value="extrapolate")
    fI = interp1d(m, Is, kind="linear", fill_value="extrapolate")
    e = imf_slope + 1.0
    sign = 1.0 if martin_sign else -1.0
    out, n_have, drawn, batch = [], 0, 0, max(2 * n_target, 1000)
    while n_have < n_target and drawn < max_factor * n_target:
        r = rs.random_sample(batch)
        ms = (r * (m_hi ** e - m_lo ** e) + m_lo ** e) ** (1 / e)
        v, i = fV(ms), fI(ms)
        if binary_fraction > 0:
            bb = rs.random_sample(batch) <= binary_fraction
            m2 = (0.5 + 0.5 * rs.random_sample(int(bb.sum()))) * ms[bb]
            v[bb] = -2.5 * np.log10(10 ** (-0.4 * v[bb]) + 10 ** (-0.4 * fV(m2)))
            i[bb] = -2.5 * np.log10(10 ** (-0.4 * i[bb]) + 10 ** (-0.4 * fI(m2)))
        rad = ast.radii(rs, batch)
        okv, ev = ast.band["F606W"].draw(rs, rad, v - ast.zp["F606W"])
        oki, ei = ast.band["F814W"].draw(rs, rad, i - ast.zp["F814W"])
        v, i = v + sign * ev, i + sign * ei
        keep = okv & oki & (v > lo) & (v < hi)
        if cut == "isochrone" and len(ri) > 1:
            v_iso, i_iso = lininterp(ri, rv, i), lininterp(rv, ri, v)
            keep &= (v <= v_iso + width) & (v >= v_iso - width) & (i <= i_iso + width) & (i >= i_iso - width)
        elif cut == "fiducial" and fid is not None:
            keep &= fid.keep(v, i)
        elif cut == "ridge":
            keep &= ridge.distance(v - i, v) <= width
        out.append(np.column_stack([v[keep], i[keep]]))
        n_have += int(keep.sum())
        drawn += batch
        rate = max(n_have, 1) / drawn
        batch = int(min(max(1000, 1.2 * (n_target - n_have) / rate), 2e6))
    cat = np.vstack(out)
    return cat[:n_target]
