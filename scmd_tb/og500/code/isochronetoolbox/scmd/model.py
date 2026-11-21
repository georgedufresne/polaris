"""Synthetic CMDs from one isochrone set: photometry, quota sampling, scatter.

A condensation of the stage-1 script of the old pipeline
(``scmd_grid_red_dist_chaboyer_g_err_3err_test.py``).  Every floating-point
operation and every random draw is made in the same order as there, so for the
same seed the catalogues agree with it to the bit (see ``tests/scmd``).  What
moved is where the work happens: the bolometric corrections depend only on the
age, so they are evaluated once per age instead of once per (DM, E(B-V)) cell,
and nothing is written to disk.

The dead code of the original (1D photometric-error CDFs, perpendicular RP-bin
errors, the single-region scatter, ``simulate_stars`` and its zero binary
fraction) never ran in the pipeline and is not carried over.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator, interp1d

ISO_COLUMNS = [
    "EEP", "isochrone_age_yr", "initial_mass", "star_mass", "star_mdot",
    "he_core_mass", "c_core_mass", "log_L", "log_LH", "log_LHe", "log_Teff",
    "log_R", "log_g", "surface_h1", "surface_he3", "surface_he4",
    "surface_c12", "surface_o16", "log_center_T", "log_center_Rho",
    "center_gamma", "center_h1", "center_he4", "center_c12",
]
BANDS = ("G", "BP", "RP")
A_PER_EBV = {"G": 2.740, "BP": 3.374, "RP": 2.035}   # A_band / E(B-V), Gaia Collab. 2018

# Normalisation of the perpendicular-scatter clip: the fixed extent of the
# NGC 6791 reference CMD, as hard-coded in the original.
COLOR_SCALE = 1.7554049999999997 - 0.9020960000000002
MAG_SCALE = 18.999708 - 15.588097


def read_isochrones(path: str) -> pd.DataFrame:
    """Every isochrone in a MIST-format file, read exactly as the original did.

    ``skiprows=10`` stops one line short of the column-name comment, so pandas
    takes the first *data* row as the header: the lowest-mass point of the
    youngest isochrone is dropped.  It is kept that way on purpose (it lies
    far below the fitted magnitude range).
    """
    df = pd.read_csv(path, sep=r"\s+", comment="#", skiprows=10)
    df.columns = ISO_COLUMNS
    return df


def select_age(iso: pd.DataFrame, age_gyr: float, exact: bool = False) -> Tuple[float, pd.DataFrame]:
    """The isochrone closest to ``age_gyr``: (age in yr, its rows).

    The original matched with ``np.isclose(..., rtol=0.01)``.  From 10.0 Gyr up
    a 0.1 Gyr neighbour is within 1 %, so those "isochrones" are the union of
    two or three ages.  ``exact=True`` takes only the closest age block.
    """
    ages = np.unique(iso["isochrone_age_yr"])
    age = ages[np.argmin(np.abs(ages - age_gyr * 1e9))]
    col = iso["isochrone_age_yr"]
    rows = iso[col == age] if exact else iso[np.isclose(col, age, rtol=0.01)]
    return float(age), rows


class BCTable:
    """MIST bolometric corrections for Gaia G, BP, RP on (Teff, log g).

    Linear interpolation on a Delaunay triangulation of the raw (Teff, log g)
    points, as in the original.  ``data/ngc6791/bc_mist_gaia_edr3_*.npz`` hold
    the tables (extracted with fidanka's loader, Av = 0, Rv = 3.1).
    """

    def __init__(self, path: str):
        z = np.load(path)
        self.feh = float(z["feh"])
        self._interp = LinearNDInterpolator(z["points"], z["bc"])

    def absolute(self, rows: pd.DataFrame) -> dict:
        """Absolute magnitudes M_bol - BC of the isochrone points."""
        teff = 10 ** rows["log_Teff"].values
        mbol = 4.75 - 2.5 * rows["log_L"].values
        bc = self._interp(teff, rows["log_g"].values)
        return {b: mbol - bc[:, k] for k, b in enumerate(BANDS)}


def apparent(absmag: dict, dm: float, ebv: float) -> dict:
    """Apparent magnitudes for an apparent distance modulus ``dm`` (in G).

    Reddening is added band by band and the true modulus is dm - A_G, in that
    order, so the sums round exactly as in the original.
    """
    true_dm = dm - A_PER_EBV["G"] * ebv
    return {b: (absmag[b] + A_PER_EBV[b] * ebv) + true_dm for b in BANDS}


@dataclass
class Quota:
    """How many stars to draw, and the completeness to thin them with.

    The observed stars with ``g_lo <= G <= g_hi`` are counted in the 0.5-mag
    completeness bins whose centres lie in that range and divided by the
    completeness; the sum is the number of stars drawn for every sCMD.
    """
    n: int
    bins: List[Tuple[float, float]]          # (centre, completeness)
    g_lo: float = 16.25
    g_hi: float = 18.0

    @classmethod
    def from_counts(cls, obs_g: np.ndarray, completeness: pd.DataFrame,
                    g_lo: float = 16.25, g_hi: float = 18.0) -> "Quota":
        centers = completeness["bin_center"].values
        frac = completeness["completeness"].values
        edges = np.array([c - 0.25 for c in centers] + [centers[-1] + 0.25])
        g = obs_g[(obs_g >= g_lo) & (obs_g <= g_hi)]
        corrected, bins = [], []
        for i, c in enumerate(centers):
            if c < g_lo or c > g_hi:
                continue
            n_obs = np.sum((g >= edges[i]) & (g < edges[i + 1]))
            corrected.append(n_obs / frac[i] if frac[i] > 0 else 0)
            bins.append((c, frac[i]))
        return cls(int(np.sum(corrected)), bins, g_lo, g_hi)


@dataclass
class Scatter:
    """Perpendicular scatter per G bin for one CMD region (red or blue)."""
    left: np.ndarray
    right: np.ndarray
    std: np.ndarray

    @classmethod
    def read(cls, path: str) -> "Scatter":
        d = pd.read_csv(path)
        return cls(d["bin_left"].values, d["bin_right"].values, d["perp_std"].values)


def _slope(g: np.ndarray, color: np.ndarray):
    """d(color)/dG along the isochrone, from central differences in G order."""
    o = np.argsort(g)
    g, color = g[o], color[o]
    dg = g[2:] - g[:-2]
    ok = np.abs(dg) > 1e-6
    s = (color[2:] - color[:-2])[ok] / dg[ok]
    x = g[1:-1][ok]
    if len(s):
        x, s = np.r_[g[0], x, g[-1]], np.r_[s[0], s, s[-1]]
    else:
        x, s = np.array([g[0], g[-1]]), np.zeros(2)
    return interp1d(x, s, kind="linear", bounds_error=False, fill_value="extrapolate")


def _scatter_region(rs, G, BP, RP, sc: Scatter, slope, max_norm: float):
    """Scatter perpendicular to the isochrone, clipped in normalised CMD units.

    Stars outside every bin of the region are dropped, as in the original.
    """
    k = np.full(len(G), -1, dtype=int)
    for i in range(len(sc.left)):
        k[(G >= sc.left[i]) & (G < sc.right[i])] = i
    v = k >= 0
    if not v.any():
        return None
    G, BP, RP = G[v], BP[v], RP[v]
    s = slope(G)
    perp = np.where(np.abs(s) > 1e-6, -1.0 / s, np.sign(s) * 1e6)
    norm = np.sqrt(1 + perp ** 2)
    err = rs.normal(0, sc.std[k[v]])
    dg, dc = err * (1.0 / norm), err * (perp / norm)
    dist = np.sqrt((dg / MAG_SCALE) ** 2 + (dc / COLOR_SCALE) ** 2)
    f = np.where(dist > max_norm, max_norm / dist, 1.0)
    dg, dc = dg * f, dc * f
    return G + dg, BP + dc * 0.5, RP + -dc * 0.5


def make_scmd(rs: np.random.RandomState, mass: np.ndarray, phot: dict, quota: Quota,
              red: Scatter, blue: Scatter, imf_slope: float = -1.02,
              max_norm: float = 0.05) -> Optional[np.ndarray]:
    """One synthetic CMD as an (n, 3) array of observed (G, BP, RP), or None.

    1. Draw ``quota.n`` masses from a power-law IMF between the lowest and
       highest isochrone mass with ``g_lo <= G <= g_hi``, and interpolate the
       magnitudes linearly in mass.
    2. Thin each completeness bin to round(n * completeness) stars.
    3. Scatter the red region (BP-RP >= 1.15 and G < 17.5) and the blue region
       perpendicular to the isochrone with their own per-bin widths.

    Returns None where the original failed (no isochrone point in the window).
    Draws: one uniform per star, one permutation per thinned bin, one normal
    per scattered star, red region first -- the original's order exactly.
    """
    G, BP, RP = phot["G"], phot["BP"], phot["RP"]
    o = np.argsort(mass)
    m, Gs = mass[o], G[o]
    win = (Gs >= quota.g_lo) & (Gs <= quota.g_hi)
    if not win.any():
        return None
    m_lo, m_hi = m[win].min(), m[win].max()
    f = [interp1d(m, y, kind="linear", fill_value="extrapolate") for y in (Gs, BP[o], RP[o])]

    e = imf_slope + 1.0
    r = rs.random_sample(quota.n)
    m_star = (r * (m_hi ** e - m_lo ** e) + m_lo ** e) ** (1 / e)
    g, bp, rp = (fi(m_star) for fi in f)

    keep = np.ones(len(g), dtype=bool)
    for c, comp in quota.bins:
        idx = np.where((g >= c - 0.25) & (g < c + 0.25))[0]
        if len(idx):
            n_rm = len(idx) - int(round(len(idx) * comp))
            if n_rm > 0:
                keep[rs.choice(idx, size=n_rm, replace=False)] = False
    g, bp, rp = g[keep], bp[keep], rp[keep]

    slope = _slope(G, BP - RP)
    is_red = (bp - rp >= 1.15) & (g < 17.5)
    parts = [_scatter_region(rs, g[m_], bp[m_], rp[m_], sc, slope, max_norm)
             for m_, sc in ((is_red, red), (~is_red, blue)) if m_.any()]
    parts = [p for p in parts if p is not None]
    if not parts:
        return np.empty((0, 3))
    return np.column_stack([np.concatenate(c) for c in zip(*parts)])
