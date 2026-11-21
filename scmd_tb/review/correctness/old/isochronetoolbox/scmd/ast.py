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
* the sCMD gets the observed fit stars' selection, and as many stars as the
  observed sample has.  ``cut="ridge"``: the observed stars are those within a
  normalised distance of a ridge line traced through the raw photometry
  (``ridge.py``), and each sCMD is cut the same way around a ridge traced
  through its own stars.  ``cut="isochrone"``: M. Ying's 0.08-mag cut in both
  bands around the sCMD's own isochrone (the observed stars were cut around
  his fiducial line).  Cutting the sCMD around the *observed* fiducial
  (``cut="fiducial"``) is a trap with a fixed star count: a misplaced
  isochrone is filled up with its error-tail stars.

The statistics afterwards are the NGC 6791 ones (``stats.py``).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline, interp1d

from .ridge import RidgeCut

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


class NearestErrors:
    """Completeness and errors from the artificial stars nearest to each synthetic star.

    For a synthetic star with true (F606W, F814W) and distance from the centre:
    the ``k`` artificial stars nearest in (input F606W, input F814W, radius),
    with magnitudes divided by ``scale_mag`` and radius by ``scale_r``; one of
    them at random; the synthetic star is lost if that artificial star was not
    recovered (0.5 px, 0.75 mag), otherwise it takes that star's own errors in
    both bands.  For speed the neighbourhood is that of the artificial star
    nearest to the synthetic star (each artificial star's ``k`` neighbours are
    found once), at most one artificial-star spacing off.  The two bands' errors of one artificial star are strongly
    correlated (a blend brightens both; r = 0.84 for 47 Tuc), which per-band
    draws (MakeTestCMD) lose, roughly doubling the colour scatter.
    """

    def __init__(self, a: pd.DataFrame, guess=(3000.0, 3000.0), guess_radius=200.0,
                 k: int = 30, scale_mag: float = 0.1, scale_r: float = 200.0):
        from scipy.spatial import cKDTree
        x, y = a["InputX"].values, a["InputY"].values
        vin, iin = a["InputF606W"].values.astype(float), a["InputF814W"].values.astype(float)
        vout, iout = a["OutputF606W"].values.astype(float), a["OutputF814W"].values.astype(float)
        xc, yc = centre(x, y, guess, guess_radius)
        r = np.hypot(x - xc, y - yc)
        self.rec = ((np.hypot(x - a["OutputX"].values, y - a["OutputY"].values) < 0.5)
                    & (np.abs(vin - vout) < 0.75) & (np.abs(iin - iout) < 0.75))
        self.ev, self.ei = vin - vout, iin - iout              # (input - output), as MakeTestCMD
        self.k, self.scale = k, np.array([scale_mag, scale_mag, scale_r])
        pts = np.column_stack([vin, iin, r]) / self.scale
        self.tree = cKDTree(pts)
        self.nbr = self.tree.query(pts, k=k)[1].astype(np.int32)

    def draw(self, rs: np.random.RandomState, radius, v_inst, i_inst):
        """(recovered, F606W error, F814W error), errors as (input - output)."""
        _, near = self.tree.query(np.column_stack([v_inst, i_inst, radius]) / self.scale)
        pick = self.nbr[near, np.minimum((rs.random_sample(len(near)) * self.k).astype(int), self.k - 1)]
        return self.rec[pick], self.ev[pick], self.ei[pick]


class ArtificialStars:
    """Both bands' error model plus the radial distribution of the observed stars."""

    def __init__(self, ast_csv: str, fit_stars: pd.DataFrame, zp: Dict[str, float],
                 window: Dict[str, Tuple[float, float]], short_long: Dict[str, float],
                 n_rad=10, n_mag=12, guess=(3000.0, 3000.0), guess_radius=200.0,
                 rows: Optional[np.ndarray] = None):
        """``rows``: a boolean mask over the artificial stars to use (default all)."""
        a = pd.read_csv(ast_csv, dtype=F32)
        if rows is not None:
            a = a[np.asarray(rows, bool)].reset_index(drop=True)
        self.n_ast = len(a)
        self._args = (ast_csv, fit_stars, zp, window, short_long, n_rad, n_mag, guess, guess_radius)
        self.zp = zp
        self.band = {}
        for b, (ci, co) in {"F606W": ("InputF606W", "OutputF606W"),
                            "F814W": ("InputF814W", "OutputF814W")}.items():
            lo, hi = window[b]
            self.band[b] = bin_artificial_stars(
                a["InputX"].values, a["InputY"].values, a[ci].values,
                a["OutputX"].values, a["OutputY"].values, a[co].values,
                lo - zp[b], hi - zp[b], short_long[b] - zp[b], n_rad, n_mag, guess, guess_radius)
        self.nearest = NearestErrors(a, guess, guess_radius)
        # RadialDensity.f90: distances of the observed fit stars from their own centre
        x, y = fit_stars["x"].values, fit_stars["y"].values
        xc, yc = centre(x, y, guess, guess_radius)
        self.dist = np.sort(np.sqrt((x.astype(F32) - F32(xc)) ** 2 + (y.astype(F32) - F32(yc)) ** 2)).astype(float)
        self.cdf = np.arange(1, len(self.dist) + 1) / len(self.dist)

    def halves(self, seed: int = 0) -> Tuple["ArtificialStars", "ArtificialStars"]:
        """The same model built from two disjoint random halves of the artificial stars,
        so that two sCMDs can be given independent photometric errors (bootstrap.py)."""
        a = self._args
        n = len(pd.read_csv(a[0], usecols=["InputX"]))
        half = np.random.RandomState(seed).random_sample(n) < 0.5
        return ArtificialStars(*a, rows=half), ArtificialStars(*a, rows=~half)

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
                  rcut: Optional[RidgeCut] = None, errors: str = "independent") -> Optional[np.ndarray]:
    """One synthetic CMD with artificial-star errors: an (n, 2) array of observed
    (F606W, F814W), n = ``n_target`` unless the cell is hopeless; None if the
    isochrone misses the window.

    Masses come from the power-law mass function between the lowest and highest
    isochrone mass with F606W within ``pad`` mag of the window, so stars can
    scatter in from outside it.  Each star gets a distance from the centre and,
    in each band, a completeness test and an artificial-star error
    (``errors="independent"``, MakeTestCMD), or the recovery and both errors of
    one of the artificial stars nearest in input magnitudes and radius
    (``"nearest"``, ``NearestErrors``).  Then the
    selection of the observed fit stars, the F606W window and

    * ``cut="ridge"``: a ridge traced, with the observed ridge's settings,
      through the first ``rcut.n_trace`` stars in ``rcut``'s trace box (as many
      as the observed trace box held), and the stars in the window and colour
      box within ``rcut.width`` of it;
    * ``cut="isochrone"``: 0.08 mag in both bands around the sCMD's own
      isochrone, as in MakeTestCMD;
    * ``cut="fiducial"``: around the observed fiducial.

    Stars are drawn in batches until ``n_target`` survive, up to ``max_factor`` x
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
    fV = interp1d(m, Vs, kind="linear", fill_value="extrapolate")
    fI = interp1d(m, Is, kind="linear", fill_value="extrapolate")
    e = imf_slope + 1.0
    sign = 1.0 if martin_sign else -1.0

    def draw(n):
        """n stars: observed (F606W, F814W) and whether they were recovered inside the window."""
        r = rs.random_sample(n)
        ms = (r * (m_hi ** e - m_lo ** e) + m_lo ** e) ** (1 / e)
        v, i = fV(ms), fI(ms)
        if binary_fraction > 0:
            bb = rs.random_sample(n) <= binary_fraction
            m2 = (0.5 + 0.5 * rs.random_sample(int(bb.sum()))) * ms[bb]
            v[bb] = -2.5 * np.log10(10 ** (-0.4 * v[bb]) + 10 ** (-0.4 * fV(m2)))
            i[bb] = -2.5 * np.log10(10 ** (-0.4 * i[bb]) + 10 ** (-0.4 * fI(m2)))
        rad = ast.radii(rs, n)
        if errors == "nearest":
            ok, ev, ei = ast.nearest.draw(rs, rad, v - ast.zp["F606W"], i - ast.zp["F814W"])
            return v + sign * ev, i + sign * ei, ok
        okv, ev = ast.band["F606W"].draw(rs, rad, v - ast.zp["F606W"])
        oki, ei = ast.band["F814W"].draw(rs, rad, i - ast.zp["F814W"])
        v, i = v + sign * ev, i + sign * ei
        return v, i, okv & oki

    def fill(select, n_need):
        """Draw in batches until ``n_need`` recovered stars pass ``select``; them, in order drawn."""
        out, n_have, drawn, batch = [], 0, 0, max(2 * n_need, 1000)
        while n_have < n_need and drawn < max_factor * n_need:
            v, i, keep = draw(batch)
            keep &= select(v, i)
            out.append(np.column_stack([v[keep], i[keep]]))
            n_have += int(keep.sum())
            drawn += batch
            rate = max(n_have, 1) / drawn
            batch = int(min(max(1000, 1.2 * (n_need - n_have) / rate), 2e6))
        return np.vstack(out)

    if cut == "ridge":
        pool = fill(lambda v, i: rcut.in_trace_box(v - i, v), rcut.n_trace)
        if len(pool) < 1000:                         # hopeless cell: too few stars for a ridge
            return pool[rcut.in_box(pool[:, 0] - pool[:, 1], pool[:, 0])][:n_target]
        t = pool[:rcut.n_trace]                      # as many as the observed trace box held
        ridge = rcut.trace(t[:, 0] - t[:, 1], t[:, 0])

        def select(v, i):
            return rcut.select(v - i, v, ridge)
        cat = pool[select(pool[:, 0], pool[:, 1])]
        if len(cat) < n_target:
            cat = np.vstack([cat, fill(select, n_target - len(cat))])
        return cat[:n_target]

    def window(v):
        return (v > lo) & (v < hi)
    if cut == "isochrone" and len(ri) > 1:
        def select(v, i):
            v_iso, i_iso = lininterp(ri, rv, i), lininterp(rv, ri, v)
            return window(v) & ((v <= v_iso + width) & (v >= v_iso - width) & (i <= i_iso + width) & (i >= i_iso - width))
    elif cut == "fiducial" and fid is not None:
        def select(v, i):
            return window(v) & fid.keep(v, i)
    else:
        def select(v, i):
            return window(v)
    return fill(select, n_target)[:n_target]
