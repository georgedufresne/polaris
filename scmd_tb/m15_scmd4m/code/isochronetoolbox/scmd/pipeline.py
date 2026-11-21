"""The sCMD grid: generate, score and keep the best, one (run, age) at a time.

The old pipeline ran three scripts per isochrone set: stage 1 wrote 308 CSVs
per age (a 28 x 11 grid in apparent distance modulus and E(B-V)), stage 2
scored each against the observed CMD with a KDE Jensen-Shannon distance and
deleted all but the best four, and stage 3 computed the 2D KS statistic (and a
second JS distance) for the survivors.  Here the three stages are one loop in
memory, with the same selection rule, so a (run, age) costs ~35 s instead of
~5 min and leaves two small CSVs instead of 308 files.

Randomness is seeded per (run, age), so any sCMD can be regenerated exactly
(:func:`regenerate`) without having been stored.
"""

from __future__ import annotations

import glob
import os
import shutil
import socket
import time
import traceback
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .model import (BCTable, Quota, Scatter, apparent, make_scmd, read_dsep, read_isochrones,
                    select_age)
from .stats import KDE, js_adaptive, js_fixed, ks_2d, ks_2d_count, normalise, rank2d

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data", "ngc6791")
DATA_TUC = os.path.join(ROOT, "data", "47tuc")
MARTIN = "/dartfs-hpc/scratch/f007058/Files_from_Martin"
A_HST = {"F606W": 2.88, "F814W": 1.86}   # A/E(B-V), ACS/WFC, MIST at the 47 Tuc turnoff, Rv = 3.1

KS_COLUMNS = ["ks_statistic", "js_divergence", "mean_diff", "p95_diff", "n_observed",
              "n_simulated", "simulation_id", "filename", "age_gyr", "distance_modulus",
              "extinction_ebv", "file_path", "plot_path", "kde_js", "kde_rank"]
JS_COLUMNS = ["simulation_id", "age_gyr", "distance_modulus", "extinction_ebv",
              "js_distance", "n_stars"]


@dataclass
class Setup:
    """Everything that defines a grid.  The defaults are the NGC 6791 pipeline
    as it ran (``test_scmd_chaboyer_account_physics.batch``); ``PRESETS`` holds
    the others (``preset("47tuc")``)."""
    name: str = "NGC6791"
    model: str = "gaia"                   # "gaia": regional scatter; "ast": artificial stars
    iso: str = "/dartfs-hpc/rc/lab/C/ChaboyerB/Isochrones_Only_July_10/{run}/isos/isochrones_example.txt"
    iso_format: str = "mist"              # "mist": MIST/pyeep + BCs; "dsep": M. Ying's DSEP files
    obs: str = os.path.join(DATA, "obs_referee_chaboyer_DM13.280_EBV0.200_age8.3Gyr_001.csv")
    obs_color: str = "BP_RP"
    obs_mag: str = "Gmag"
    completeness: str = os.path.join(DATA, "completeness.csv")
    scatter_red: str = os.path.join(DATA, "scatter_red.csv")
    scatter_blue: str = os.path.join(DATA, "scatter_blue.csv")
    bc: str = os.path.join(DATA, "bc_mist_gaia_edr3_fehm300.npz")
    ages: Tuple[float, ...] = tuple(round(4.0 + 0.1 * i, 1) for i in range(81))
    dm_grid: Tuple[float, float, float] = (13.23, 13.51, 0.01)     # np.arange(*dm_grid)
    ebv_grid: Tuple[float, float, float] = (0.10, 0.21, 0.01)
    exact_age: bool = False               # True: no 1 % age union above 10 Gyr
    keep: int = 4                         # sCMDs per age passed from stage 2 to stage 3
    chunk: int = 9                        # ages per unit of work (part of the unit names)
    kde_window: Tuple[Tuple[float, float], Tuple[float, float]] = ((1.0, 1.5), (15.0, 19.0))
    ks_window: Tuple[Tuple[float, float], Tuple[float, float]] = ((0.9, 1.7), (16.25, 18.0))
    ks_3sigma: float = 0.0468             # bootstrap null: mean 0.0268 + 3 x 0.0067
    ks_model_ecdf: str = "interp"         # stage-3 model ECDF at the observed stars: interp (linear interpolation
                                          # of the sCMD ECDF, production) | count (exact count; large n_model)
    ks_percentile: float = 95.0           # stage-3 2D KS = this percentile of |ECDF diff| (p95_diff keeps the 95th)
    imf_slope: float = -1.02
    # artificial-star model (model = "ast"); see ast.py
    ast: str = ""                         # artificial-star catalogue
    fit_stars: str = ""                   # observed fit stars: radial distribution, AST windows
    fiducial: str = ""                    # observed fiducial line (v, vi)
    zp_v: float = 0.0                     # VEGAmag - instrumental, F606W
    zp_i: float = 0.0                     # VEGAmag - instrumental, F814W
    short_long_v: float = 0.0             # short/long exposure boundary, F606W (VEGAmag)
    short_long_i: float = 0.0
    sgb_v: float = 0.0                    # observed subgiant branch, F606W
    window_dmag: float = 2.0              # fit window: sgb_v +/- window_dmag
    pad: float = 1.0                      # stars are drawn this far outside the window
    fid_width: float = 0.08               # half-width of the cut around the fiducial line
    binary_fraction: float = 0.0
    martin_sign: bool = False             # True: observed = true + (input - output)
    ast_errors: str = "independent"       # independent: per-band draws (MakeTestCMD) | nearest: copy both
                                          # errors of one of the ast_k nearest artificial stars | smooth:
                                          # draw from a kernel density of those neighbours' errors (ast.py)
    ast_k: int = 30                       # artificial stars per neighbourhood (nearest, smooth)
    ast_metric: str = "linear"            # neighbour search: linear (0.1 mag, 200 px) | crowding (0.1 mag,
                                          # ln(r + 200 px) / 0.09, short/long exposures kept apart)
    ast_split: bool = True                # ast_metric="crowding": keep short- and long-exposure artificial
                                          # stars in separate neighbourhoods (hard split at short_long_v/i);
                                          # False: neighbourhoods span the 16.25-16.75 transition, so the drawn
                                          # artificial star's own exposure regime decides (47 Tuc: deep-frame
                                          # share 22% at F606W 16.25-16.5, 80% at 16.5-16.75)
    radial: str = "observed"              # radii: observed fit-star distances | corrected: per 0.25-mag
                                          # bin, weighted by 1/completeness (ast.py)
    radial_source: str = "fit"            # radial="corrected" built from: fit = the fit stars (after the cut,
                                          # so the cut's core losses are applied twice) | raw = obs_raw stars
                                          # within radial_loose mag per band of the fiducial, over the generated
                                          # range (window +/- pad)
    radial_loose: float = 0.3
    radial_cfloor: float = 0.1            # floor of the completeness in the 1/c weights
    ast_scale_mag: float = 0.1            # nearest/smooth neighbourhoods: magnitude scale
    ast_scale_lnr: float = 0.09           # crowding metric: scale of ln(r + 200 px)
    crn: bool = False                     # common random numbers: every cell of every (run, age) drawn
    crn_seed: int = 20261025              # from RandomState(crn_seed), so cells differ only by parameters
    extra_sigma_v: float = 0.0            # PSF-mismatch scatter: independent Gaussian errors added to every
    extra_sigma_i: float = 0.0            # synthetic star in F606W and F814W, on top of the AST error
    extra_sigma_table: tuple = ()         # or ((F606W, sigma), ...): sigma interpolated at each star's true
                                          # F606W, used in both bands (independent draws); overrides the two above
    zp_correct: bool = False              # cut="ridge": remove the observed colour zero-point map (local median
    zp_k: int = 50                        # colour offset of the zp_k nearest upper-MS stars, minus the offset
    zp_ref: Tuple[float, float] = (17.6, 18.9)   # the ASTs predict there) before the ridge selection
    n_model: int = 0                      # stars in each model sCMD (0: as many as the observed fit stars);
                                          # the ridge is still traced on an observed-sized sample
    zp_mag_frac: float = 0.0              # share of the colour correction taken by F606W (0: F814W only;
                                          # M15's pattern: dF606W = 0.57 dcolour, diag/diag3.py)
    cut: str = "isochrone"                # ridge | isochrone (0.08 mag) | fiducial | none
    obs_raw: str = ""                     # cut="ridge": the raw photometry the fit stars come from
    ridge_width: float = 0.05             # cut="ridge": normalised perpendicular distance (ridge.py)
    ridge_box_color: Tuple[float, float] = (0.35, 0.95)   # colour box the ridge is traced in
    ridge_margin: float = 0.3             # ridge traced over the window +/- this, then trimmed to it
    ridge_exclude: tuple = ()             # ((colour range), (mag range)) boxes left out, observed and sCMD
    ridge_scales: Optional[Tuple[float, float]] = None    # axis normalisation; None: observed range in the box
    ridge_warp: Optional[tuple] = None    # (c0, ((colour, dmag), ...), cap), see ridge.build_smooth_warp
    ridge_warp_from_turnoff: bool = False  # warp colours are offsets from each CMD's turnoff colour
    ridge_smoothing: float = 3e-5         # smoothing spline: s = this x ridge points (raw CMD units)

    reddening: str = "ebv"                # grid in E(B-V) (band coefficients) or "evi":
                                          # E(V-I) colour shift, DM applied to F606W (M. Ying)

    def __post_init__(self):
        self.dm = np.arange(*self.dm_grid)
        self.ebv = np.arange(*self.ebv_grid)
        for name, allowed in (("model", ("gaia", "ast")), ("iso_format", ("mist", "dsep")),
                              ("ast_errors", ("independent", "nearest", "smooth")),
                              ("ast_metric", ("linear", "crowding")), ("radial", ("observed", "corrected")), ("radial_source", ("fit", "raw")),
                              ("cut", ("ridge", "isochrone", "fiducial", "none")), ("reddening", ("ebv", "evi"))):
            if getattr(self, name) not in allowed:
                raise ValueError(f"Setup.{name} = {getattr(self, name)!r}; expected one of {allowed}")


PRESETS = {
    "ngc6791": {},
    # 47 Tuc, HST ACS F606W/F814W (M. Ying's data).  mu 13.15-13.45 (Denissenkov et
    # al. 2017, Brogaard et al. 2017) applied to F606W, E(V-I) 0-0.05 (Denissenkov
    # et al. 2017) as a colour shift -- the ranges and convention of his global_var.py
    "47tuc": dict(
        name="47Tuc", model="ast",
        iso="/dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso/{run}/isos/isochrones_example.txt",
        obs=f"{MARTIN}/47Tuc_data/47Tuc_fitstars_2dKS.dat", obs_color="vi", obs_mag="v",
        fit_stars=f"{MARTIN}/47Tuc_data/47Tuc_fitstars.dat",
        ast=f"{MARTIN}/47Tuc_data/47Tucartstars.dat",
        fiducial=f"{MARTIN}/47Tuc_data/fiducial_lines.csv",
        bc=os.path.join(DATA_TUC, "bc_mist_acs_wfc_fehm075.npz"),
        zp_v=30.5742, zp_i=29.6777, short_long_v=16.48, short_long_i=15.55, sgb_v=17.02,
        ages=tuple(round(8.0 + 0.1 * i, 1) for i in range(61)),   # the tracks support <= 14 Gyr
        dm_grid=(13.15, 13.455, 0.01), ebv_grid=(0.0, 0.0525, 0.005), reddening="evi",
        kde_window=((0.3, 1.3), (15.02, 19.02)), ks_window=((0.3, 1.3), (15.02, 19.02)),
        # present-day mass function dN/dm ~ m^alpha, alpha = -0.45 +- 0.12 (Ebrahimi et al. 2020,
        # Table 2: global slope over 0.2-0.8 Msun; M. Ying's global_var.py has -0.84)
        imf_slope=-0.45, exact_age=True,
        # errors: drawn from a kernel density of the 30 artificial stars nearest in (V_in, I_in,
        # radius), searched in crowding-scaled units with short and long exposures apart; radii from
        # completeness-corrected distributions per magnitude; the same random numbers in every cell;
        # binaries 2 % (q = 0.5-1)
        ast_errors="smooth", ast_k=30, ast_metric="crowding", radial="corrected", crn=True,
        binary_fraction=0.02,
        # observed colours corrected for the zero-point map (rms 0.007); PSF-mismatch scatter added to
        # the sCMDs, sigma per band by F606W, calibrated so that the model main sequence has the
        # corrected data's core width (2026-10-25)
        zp_correct=True,
        extra_sigma_table=((16.0, 0.0074), (17.8, 0.0031), (18.2, 0.004), (18.6, 0.0049), (18.9, 0.0059)),
        keep=8,                           # stage 2 (KDE, bandwidth 0.15) passes 8 cells per age to the KS
        ks_3sigma=0.0105,                 # 47 Tuc bootstrap 3 sigma (fresh references, resampled ASTs);
                                          # only collect's filtered table uses it
        # own fit stars (ridge.py): the raw photometry in the window and colour box, a
        # ridge traced through it with the NGC 6791 method, stars within 0.15 of it in
        # normalised perpendicular distance (NGC 6791 used 0.05); every sCMD likewise.
        # Warp: the NGC 6791 knots with their colour offsets halved (47 Tuc's subgiant
        # branch is ~0.1 wide in colour, NGC 6791's ~0.2), placed from each CMD's turnoff
        cut="ridge", obs_raw=f"{MARTIN}/47Tuc_data/47Tuc.dat", ridge_width=0.15,
        ridge_box_color=(0.35, 0.95), ridge_exclude=(((0.40, 0.60), (15.5, 16.75)),),   # blue stragglers / noise
        ridge_warp=(0.04, ((0.055, -0.05), (0.095, -0.15), (0.135, -0.30), (0.19, -0.45)), -0.5),
        ridge_warp_from_turnoff=True, ridge_smoothing=5e-5),
}


def preset(name: str, **override) -> Setup:
    return Setup(**{**PRESETS[name.lower()], **override})


def _window(pts: np.ndarray, w) -> np.ndarray:
    (c0, c1), (m0, m1) = w
    return (pts[:, 0] >= c0) & (pts[:, 0] <= c1) & (pts[:, 1] >= m0) & (pts[:, 1] <= m1)


def seed(run: int, age: float) -> int:
    return int(np.random.SeedSequence([int(run), int(round(age * 10))]).generate_state(1)[0])


class Context:
    """The per-process state: inputs read once, KDE grid and observed pdf built once."""

    def __init__(self, s: Setup):
        self.s = s
        self.bc = BCTable(s.bc)
        obs = pd.read_csv(s.obs)
        if "BP_RP" not in obs.columns and "BPmag" in obs.columns:
            obs["BP_RP"] = obs["BPmag"] - obs["RPmag"]
        obs = self._model_inputs(s, obs)
        pts = obs[[s.obs_color, s.obs_mag]].to_numpy()
        self.obs = pts
        o2 = pts[_window(pts, s.kde_window)]
        o2 = o2[np.isfinite(o2).all(axis=1)]
        self.kde = KDE(*s.kde_window)
        self.p_obs = normalise(self.kde.pdf(o2)) if len(o2) >= 10 else None
        self.obs3 = pts[_window(pts, s.ks_window)]
        self.obs3_cdf = rank2d(self.obs3) / len(self.obs3)
        self._iso: Dict[int, pd.DataFrame] = {}

    def _model_inputs(self, s: Setup, obs: pd.DataFrame) -> pd.DataFrame:
        self.quota = Quota.from_counts(obs["Gmag"].values, pd.read_csv(s.completeness))
        self.red, self.blue = Scatter.read(s.scatter_red), Scatter.read(s.scatter_blue)
        return obs


    def isochrones(self, run: int) -> pd.DataFrame:
        if run not in self._iso:
            read = read_dsep if self.s.iso_format == "dsep" else read_isochrones
            self._iso = {run: read(self.s.iso.format(run=run))}
        return self._iso[run]

    def points(self, cat: np.ndarray) -> np.ndarray:
        """(colour, magnitude) of a catalogue: (BP-RP, G) from (G, BP, RP)."""
        return np.column_stack([cat[:, 1] - cat[:, 2], cat[:, 0]])

    def model(self, run: int, age: float):
        """(age in yr, make): ``make(rs, dm, ebv)`` draws one catalogue of this
        (run, age) from the random state ``rs``."""
        age_yr, rows = select_age(self.isochrones(run), age, self.s.exact_age)
        absmag = self.bc.absolute(rows)
        mass = rows["initial_mass"].values

        def make(rs, dm, ebv):
            return make_scmd(rs, mass, apparent(absmag, dm, ebv), self.quota, self.red, self.blue,
                             self.s.imf_slope)
        return age_yr, make

    def grid(self, run: int, age: float):
        """(age in yr, iterator of (dm, ebv, catalogue or None)) in the original order: all
        cells drawn in turn from the (run, age) seed, or with ``crn`` each from a fresh
        RandomState(crn_seed) -- the same random numbers for every cell of every run and age."""
        age_yr, make = self.model(run, age)
        if self.s.crn:
            def cells():
                for dm in self.s.dm:
                    for ebv in self.s.ebv:
                        yield dm, ebv, make(np.random.RandomState(self.s.crn_seed), dm, ebv)
            return age_yr, cells()
        rs = np.random.RandomState(seed(run, age))

        def cells():
            for dm in self.s.dm:
                for ebv in self.s.ebv:
                    yield dm, ebv, make(rs, dm, ebv)
        return age_yr, cells()

    def stage2(self, cat: Optional[np.ndarray]):
        """(JS distance, n) on the fixed grid, or None where the original skipped."""
        if cat is None:
            return None
        pts = self.points(cat)
        pts = pts[_window(pts, self.s.kde_window)]
        if len(pts) < 100:
            return None
        js = js_fixed(self.kde, self.p_obs, pts) if self.p_obs is not None else np.inf
        return js, len(pts)

    def stage3(self, cat: np.ndarray) -> dict:
        """2D KS and adaptive JS of a catalogue against the observed stars."""
        return self.compare(self.obs3, self.obs3_cdf, cat)

    def compare(self, ref: np.ndarray, ref_cdf: np.ndarray, cat: np.ndarray) -> dict:
        """Stage 3 against any reference (colour, mag) set already cut to the KS
        window, with its ``rank2d / n``: the observed stars, or a reference sCMD
        in the bootstrap."""
        sim = self.points(cat)
        sim = sim[_window(sim, self.s.ks_window)]
        f = ks_2d_count if self.s.ks_model_ecdf == "count" else ks_2d
        ks, d = f(ref, sim, ref_cdf, self.s.ks_percentile)
        return {"ks_statistic": ks, "js_divergence": js_adaptive(ref, sim),
                "mean_diff": np.mean(d), "p95_diff": np.percentile(d, 95),
                "n_observed": len(ref), "n_simulated": len(sim)}



class ASTContext(Context):
    """HST clusters: artificial-star completeness and errors, sampled with the
    star's distance from the cluster centre (ast.py), in F606W / F814W."""

    def _model_inputs(self, s: Setup, obs: pd.DataFrame) -> pd.DataFrame:
        from .ast import ArtificialStars, Fiducial
        from .ridge import RidgeCut, Tracer
        lo, hi = s.sgb_v - s.window_dmag, s.sgb_v + s.window_dmag
        self.rcut = None
        if s.cut == "ridge":
            # our own fit stars: the raw photometry cut around a ridge line traced
            # through it; they also give the radial distribution
            raw = pd.read_csv(s.obs_raw, sep=r"\s+")
            raw = pd.DataFrame({"v": raw["Vvega"], "i": raw["Ivega"], "vi": raw["VIvega"],
                                "x": raw["x"], "y": raw["y"]})
            tracer = Tracer(warp=s.ridge_warp, warp_from_turnoff=s.ridge_warp_from_turnoff,
                            smoothing=s.ridge_smoothing)
            args = ((lo, hi), s.ridge_box_color, s.ridge_width, s.ridge_margin, s.ridge_scales, s.ridge_exclude)
            self.rcut, keep = RidgeCut.from_observed(tracer, raw.vi.values, raw.v.values, *args)
            self.zp_map = None
            if s.zp_correct:
                from .ast import colour_zero_point
                corr = colour_zero_point(raw, keep, self.rcut.ridge, s.ast, s.zp_v, s.zp_i, s.zp_ref, s.zp_k)
                raw["v"] = raw["v"] - s.zp_mag_frac * corr      # F606W takes zp_mag_frac of it
                raw["vi"] = raw["vi"] - corr
                raw["i"] = raw["v"] - raw["vi"]
                self.zp_map = corr
                self.rcut, keep = RidgeCut.from_observed(tracer, raw.vi.values, raw.v.values, *args)
            obs = raw[keep].reset_index(drop=True)
            fit = obs
        else:
            fit = pd.read_csv(s.fit_stars)
        inwin = fit[(fit["v"] >= lo) & (fit["v"] <= hi)]
        # AST magnitude ranges: the fit stars' extent inside the window (M. Ying's notebook)
        windows = {"F606W": (inwin["v"].min(), inwin["v"].max()), "F814W": (inwin["i"].min(), inwin["i"].max())}
        zp = {"F606W": s.zp_v, "F814W": s.zp_i}
        self.radial_stars = None
        if s.radial_source == "raw":
            # radial distribution from the photometry before the cut: a loose per-band band around the
            # fiducial (drops HB, BSS, field) over everything the generator draws (window +/- pad)
            r0 = pd.read_csv(s.obs_raw, sep=r"\s+")
            fl = pd.read_csv(s.fiducial)
            fv, fi = fl["v"].values, (fl["v"] - fl["vi"]).values
            ov, oi = np.argsort(fv), np.argsort(fi)
            v0, i0 = r0["Vvega"].values, r0["Ivega"].values
            v_need, i_need = np.interp(i0, fi[oi], fv[oi]), np.interp(v0, fv[ov], fi[ov])
            m = ((np.abs(v0 - v_need) < s.radial_loose) & (np.abs(i0 - i_need) < s.radial_loose)
                 & (v0 >= lo - s.pad) & (v0 <= hi + s.pad))
            for (c_lo, c_hi), (m_lo, m_hi) in (s.ridge_exclude or ()):   # the same boxes as the ridge cut (BSS)
                c0 = v0 - i0
                m &= ~((c0 >= c_lo) & (c0 <= c_hi) & (v0 >= m_lo) & (v0 <= m_hi))
            self.radial_stars = pd.DataFrame({"v": v0[m], "i": i0[m], "x": r0["x"].values[m], "y": r0["y"].values[m]})
        self.ast = ArtificialStars(s.ast, fit, zp, windows, {"F606W": s.short_long_v, "F814W": s.short_long_i},
                                   k=s.ast_k, smooth=s.ast_errors == "smooth", metric=s.ast_metric,
                                   radial=s.radial, split_exposures=s.ast_split,
                                   radial_stars=self.radial_stars, cfloor=s.radial_cfloor,
                                   scale_mag=s.ast_scale_mag, scale_lnr=s.ast_scale_lnr)
        self.fid = Fiducial(s.fiducial, s.fid_width) if s.fiducial else None
        self.v_window = (lo, hi)
        o = obs[(obs[s.obs_mag] > lo) & (obs[s.obs_mag] < hi)]
        self.n_target = len(o)
        return obs

    def points(self, cat: np.ndarray) -> np.ndarray:
        """(F606W-F814W, F606W) from (F606W, F814W)."""
        return np.column_stack([cat[:, 0] - cat[:, 1], cat[:, 0]])

    def absolute(self, rows: pd.DataFrame) -> dict:
        if self.s.iso_format == "dsep":
            return {"F606W": rows["F606W"].values, "F814W": rows["F814W"].values}
        bc = self.bc.absolute(rows)
        return {"F606W": bc["F606W"], "F814W": bc["F814W"]}

    def model(self, run: int, age: float):
        from .ast import make_scmd_ast
        age_yr, rows = select_age(self.isochrones(run), age, self.s.exact_age)
        absmag = self.absolute(rows)
        mass = rows["initial_mass"].values
        s = self.s
        coef = {"F606W": 1.0, "F814W": 0.0} if s.reddening == "evi" else A_HST

        def make(rs, dm, ebv, with_radius=False, n=None):
            n = n or s.n_model or self.n_target      # n_model: model sCMD size; n=self.n_target: data-sized
            return make_scmd_ast(rs, mass, apparent(absmag, dm, ebv, coef, "F606W"), self.ast, self.fid,
                                 n, self.v_window, s.pad, s.imf_slope, s.binary_fraction,
                                 s.martin_sign, s.cut, s.fid_width, rcut=self.rcut, errors=s.ast_errors,
                                 extra_sigma=(s.extra_sigma_v, s.extra_sigma_i),
                                 extra_sigma_table=s.extra_sigma_table or None, with_radius=with_radius)
        return age_yr, make


def context(s: Setup) -> Context:
    return ASTContext(s) if s.model == "ast" else Context(s)


def label(age_yr: float, dm: float, ebv: float, prefix: str = "NGC6791") -> Tuple[str, float, float, float]:
    """Legacy file name and the (age, DM, E(B-V)) parsed back from it."""
    name = f"{prefix}_age{age_yr / 1e9:.1f}Gyr_DM{dm:.3f}_EBV{ebv:.3f}.csv"
    return name, float(f"{age_yr / 1e9:.1f}"), float(f"{dm:.3f}"), float(f"{ebv:.3f}")


def process_age(ctx: Context, run: int, age: float, catalogs: Optional[dict] = None):
    """One (run, age): all grid cells through stage 2, the best ``keep`` through stage 3."""
    age_yr, cells = ctx.grid(run, age)
    js_rows, scored, written = [], [], []
    for dm, ebv, cat in cells:
        if cat is None:                       # stage 1 failed: no file was written
            continue
        name, a, d, e = label(age_yr, dm, ebv, ctx.s.name)
        written.append((np.nan, 0, name, a, d, e, cat))
        r = ctx.stage2(cat)
        if r is None:
            continue
        js_rows.append({"simulation_id": run, "age_gyr": a, "distance_modulus": d,
                        "extinction_ebv": e, "js_distance": r[0], "n_stars": r[1]})
        scored.append((r[0], name, a, d, e, cat))
    if scored:
        order = np.argsort([k[0] for k in scored], kind="stable")[: ctx.s.keep]
        survivors = [(scored[i][0], rank) + scored[i][1:] for rank, i in enumerate(order, 1)]
    else:
        # the old stage 2 returned early when no sCMD had 100 stars in its window
        # and deleted nothing, so stage 3 scored every file of the age
        survivors = written
    ks_rows = []
    for js, rank, name, a, d, e, cat in survivors:
        try:
            row = ctx.stage3(cat)
        except Exception:                     # the old stage 3 skipped such files
            continue
        row.update(simulation_id=run, filename=name, age_gyr=a, distance_modulus=d,
                   extinction_ebv=e, file_path="", plot_path="", kde_js=js, kde_rank=rank)
        ks_rows.append(row)
        if catalogs is not None:
            catalogs[name] = cat
    return js_rows, ks_rows


def regenerate(ctx: Context, run: int, age: float, dm: float, ebv: float) -> np.ndarray:
    """The catalogue of one grid cell, replayed from its (run, age) seed (with ``crn``
    drawn directly: every cell starts from the same seed)."""
    if ctx.s.crn:
        i, j = int(np.argmin(np.abs(ctx.s.dm - dm))), int(np.argmin(np.abs(ctx.s.ebv - ebv)))
        if abs(ctx.s.dm[i] - dm) >= 5e-4 or abs(ctx.s.ebv[j] - ebv) >= 5e-4:
            raise KeyError(f"no grid cell DM={dm} E(B-V)={ebv}")
        _, make = ctx.model(run, age)
        return make(np.random.RandomState(ctx.s.crn_seed), ctx.s.dm[i], ctx.s.ebv[j])   # the scored values
    for d, e, cat in ctx.grid(run, age)[1]:
        if abs(d - dm) < 5e-4 and abs(e - ebv) < 5e-4:
            return cat
    raise KeyError(f"no grid cell DM={dm} E(B-V)={ebv}")


# ---------------------------------------------------------------------------
# units of work, claims and collection
# ---------------------------------------------------------------------------

def units(runs: Sequence[int], ages: Sequence[float], chunk: int) -> List[Tuple[int, Tuple[float, ...]]]:
    return [(int(r), tuple(ages[i:i + chunk])) for r in runs for i in range(0, len(ages), chunk)]


def unit_key(run: int, ages: Sequence[float]) -> str:
    return f"{run}_{ages[0]:.1f}-{ages[-1]:.1f}"


def _write_csv(df: pd.DataFrame, path: str) -> None:
    tmp = f"{path}.tmp{os.getpid()}"
    df.to_csv(tmp, index=False)
    os.replace(tmp, path)


def run_unit(ctx: Context, out: str, run: int, ages: Sequence[float],
             save_catalogs: bool = False) -> None:
    d = os.path.join(out, "units", str(run))
    os.makedirs(d, exist_ok=True)
    key = unit_key(run, ages)
    js_all, ks_all, cats = [], [], ({} if save_catalogs else None)
    for age in ages:
        js, ks = process_age(ctx, run, age, cats)
        js_all += js
        ks_all += ks
    if cats:
        np.savez_compressed(os.path.join(d, key + ".npz"), **cats)
    _write_csv(pd.DataFrame(js_all, columns=JS_COLUMNS), os.path.join(d, key + ".js.csv"))
    _write_csv(pd.DataFrame(ks_all, columns=KS_COLUMNS), os.path.join(d, key + ".ks.csv"))


def worker(s: Setup, out: str, runs: Sequence[int], max_hours: float = 8.0,
           save_catalogs: bool = False, log=print) -> int:
    """Claim and process units until none are left or ``max_hours`` have passed.

    A unit is claimed by creating ``claims/<key>`` (``mkdir`` is atomic, also on
    NFS), so any number of workers -- on any account -- can share one run list
    without duplicating work.  A unit is done when its ``.ks.csv`` exists.

    The isochrone file is looked at before claiming, so a path this task cannot
    see (missing, or an expired Kerberos ticket) never strands a claim.  A unit
    that fails keeps its claim with the traceback in ``claims/<key>/error``
    (``status --release`` frees it); after three failures in a row the worker
    stops rather than work through the list on a broken node.
    """
    os.makedirs(os.path.join(out, "claims"), exist_ok=True)
    t0, n, fails, ctx, unseen = time.time(), 0, 0, None, set()
    for run, ages in units(runs, s.ages, s.chunk):
        if time.time() - t0 > max_hours * 3600:
            break
        key = unit_key(run, ages)
        if os.path.exists(os.path.join(out, "units", str(run), key + ".ks.csv")):
            continue
        iso = s.iso.format(run=run)
        if not os.path.exists(iso):
            if run not in unseen:
                log(f"run {run}: cannot see {iso}; skipped")
                unseen.add(run)
            continue
        claim = os.path.join(out, "claims", key)
        try:
            os.mkdir(claim)
        except FileExistsError:
            continue
        with open(os.path.join(claim, "owner"), "w") as fh:
            fh.write(f"{os.environ.get('SLURM_JOB_ID', '-')}_{os.environ.get('SLURM_ARRAY_TASK_ID', '-')} "
                     f"{socket.gethostname()} {time.strftime('%F %T')}\n")
        t = time.time()
        try:
            ctx = ctx or context(s)
            run_unit(ctx, out, run, ages, save_catalogs)
        except Exception:
            err = traceback.format_exc()
            log(f"{key}: FAILED\n{err}")
            with open(os.path.join(claim, "error"), "w") as fh:
                fh.write(err)
            fails += 1
            if fails >= 3:
                log("three units failed in a row; stopping this worker")
                break
            continue
        fails = 0
        n += 1
        log(f"{key}: {time.time() - t:.0f} s")
    return n


def status(s: Setup, out: str, runs: Sequence[int], release: bool = False) -> dict:
    """Units done, claimed without output yet, and untouched.

    ``release`` deletes the claims of units without output, so a new worker
    picks them up -- only do that when no worker is running.
    """
    done, claimed, todo = [], [], []
    for run, ages in units(runs, s.ages, s.chunk):
        key = unit_key(run, ages)
        claim = os.path.join(out, "claims", key)
        if os.path.exists(os.path.join(out, "units", str(run), key + ".ks.csv")):
            done.append(key)
        elif os.path.isdir(claim):
            owner = os.path.join(claim, "owner")
            who = open(owner).read().strip() if os.path.exists(owner) else "?"
            bad = "  FAILED, see error" if os.path.exists(os.path.join(claim, "error")) else ""
            claimed.append(f"{key} ({who}){bad}")
            if release:
                shutil.rmtree(claim)
        else:
            todo.append(key)
    return {"done": done, "claimed": claimed, "todo": todo}


def collect(s: Setup, out: str) -> Tuple[pd.DataFrame, dict]:
    """Per-run ``2d_ks_results_<run>.csv`` / ``best_fit_<run>.csv`` (the legacy
    layout, so the old merge and analysis scripts read them unchanged), plus
    ``combined_2d_ks_results.csv``, ``filtered_ks_scmd.csv`` (KS <= 3 sigma, cut
    in float32 as the old prefilter did) and every stage-2 score in
    ``kde_js_scores.csv.gz``.

    Only runs with every unit done are collected (the old pipeline wrote a
    run's table only after all its ages); a run whose sCMDs never reached
    stage 3 gets no table, as before.  Returns (combined table, report).
    """
    ks_dir = os.path.join(out, "2d_ks")
    os.makedirs(ks_dir, exist_ok=True)
    keys = [unit_key(0, s.ages[i:i + s.chunk])[2:] for i in range(0, len(s.ages), s.chunk)]
    allks, alljs, report = [], [], {"complete": [], "incomplete": [], "no_rows": []}
    for d in sorted(glob.glob(os.path.join(out, "units", "*"))):
        run = int(os.path.basename(d))
        paths = [os.path.join(d, f"{run}_{k}") for k in keys]
        if not all(os.path.exists(p + ".ks.csv") for p in paths):
            report["incomplete"].append(run)
            continue
        report["complete"].append(run)
        ks = [t for t in (pd.read_csv(p + ".ks.csv") for p in paths) if len(t)]
        alljs += [t for t in (pd.read_csv(p + ".js.csv") for p in paths) if len(t)]
        if not ks:
            report["no_rows"].append(run)
            continue
        ks = pd.concat(ks, ignore_index=True)
        allks.append(ks)
        ks.to_csv(os.path.join(ks_dir, f"2d_ks_results_{run}.csv"), index=False)
        best = pd.DataFrame([ks.loc[ks["ks_statistic"].idxmin()], ks.loc[ks["js_divergence"].idxmin()]])
        best["metric"] = ["best_ks", "best_js"]
        best.to_csv(os.path.join(ks_dir, f"best_fit_{run}.csv"), index=False)
    if not allks:
        return pd.DataFrame(columns=KS_COLUMNS), report
    comb = pd.concat(allks, ignore_index=True).sort_values(["simulation_id", "ks_statistic"])
    comb.to_csv(os.path.join(out, "combined_2d_ks_results.csv"), index=False)
    f = comb[["ks_statistic", "simulation_id", "age_gyr", "distance_modulus", "extinction_ebv"]].astype(
        {"ks_statistic": np.float32, "simulation_id": np.int32, "age_gyr": np.float32,
         "distance_modulus": np.float32, "extinction_ebv": np.float32})
    f[f["ks_statistic"] <= s.ks_3sigma].to_csv(os.path.join(out, "filtered_ks_scmd.csv"), index=False)
    if alljs:
        pd.concat(alljs, ignore_index=True).to_csv(os.path.join(out, "kde_js_scores.csv.gz"), index=False)
    return comb, report
