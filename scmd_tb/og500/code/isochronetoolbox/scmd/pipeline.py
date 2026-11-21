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
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from .model import BCTable, Quota, Scatter, apparent, make_scmd, read_isochrones, select_age
from .stats import KDE, js_adaptive, js_fixed, ks_2d, normalise, rank2d

DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                    "data", "ngc6791")

KS_COLUMNS = ["ks_statistic", "js_divergence", "mean_diff", "p95_diff", "n_observed",
              "n_simulated", "simulation_id", "filename", "age_gyr", "distance_modulus",
              "extinction_ebv", "file_path", "plot_path", "kde_js", "kde_rank"]
JS_COLUMNS = ["simulation_id", "age_gyr", "distance_modulus", "extinction_ebv",
              "js_distance", "n_stars"]


@dataclass
class Setup:
    """Everything that defines a grid.  The defaults are the NGC 6791 pipeline
    as it ran (``test_scmd_chaboyer_account_physics.batch``)."""
    iso: str = "/dartfs-hpc/rc/lab/C/ChaboyerB/Isochrones_Only_July_10/{run}/isos/isochrones_example.txt"
    obs: str = os.path.join(DATA, "obs_referee_chaboyer_DM13.280_EBV0.200_age8.3Gyr_001.csv")
    completeness: str = os.path.join(DATA, "completeness.csv")
    scatter_red: str = os.path.join(DATA, "scatter_red.csv")
    scatter_blue: str = os.path.join(DATA, "scatter_blue.csv")
    bc: str = os.path.join(DATA, "bc_mist_gaia_edr3_fehm300.npz")
    ages: Tuple[float, ...] = tuple(round(4.0 + 0.1 * i, 1) for i in range(81))
    dm: np.ndarray = field(default_factory=lambda: np.arange(13.23, 13.51, 0.01))
    ebv: np.ndarray = field(default_factory=lambda: np.arange(0.10, 0.21, 0.01))
    exact_age: bool = False               # True: no 1 % age union above 10 Gyr
    keep: int = 4                         # sCMDs per age passed from stage 2 to stage 3
    kde_window: Tuple[Tuple[float, float], Tuple[float, float]] = ((1.0, 1.5), (15.0, 19.0))
    ks_window: Tuple[Tuple[float, float], Tuple[float, float]] = ((0.9, 1.7), (16.25, 18.0))
    ks_3sigma: float = 0.0468             # bootstrap null: mean 0.0268 + 3 x 0.0067


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
        if "BP_RP" not in obs.columns:
            obs["BP_RP"] = obs["BPmag"] - obs["RPmag"]
        self.quota = Quota.from_counts(obs["Gmag"].values, pd.read_csv(s.completeness))
        self.red, self.blue = Scatter.read(s.scatter_red), Scatter.read(s.scatter_blue)
        pts = obs[["BP_RP", "Gmag"]].to_numpy()
        o2 = pts[_window(pts, s.kde_window)]
        o2 = o2[np.isfinite(o2).all(axis=1)]
        self.kde = KDE(*s.kde_window)
        self.p_obs = normalise(self.kde.pdf(o2)) if len(o2) >= 10 else None
        self.obs3 = pts[_window(pts, s.ks_window)]
        self.obs3_cdf = rank2d(self.obs3) / len(self.obs3)
        self._iso: Dict[int, pd.DataFrame] = {}

    def isochrones(self, run: int) -> pd.DataFrame:
        if run not in self._iso:
            self._iso = {run: read_isochrones(self.s.iso.format(run=run))}
        return self._iso[run]

    def grid(self, run: int, age: float):
        """(age in yr, iterator of (dm, ebv, catalogue or None)) in the original order."""
        age_yr, rows = select_age(self.isochrones(run), age, self.s.exact_age)
        absmag = self.bc.absolute(rows)
        mass = rows["initial_mass"].values
        rs = np.random.RandomState(seed(run, age))

        def cells():
            for dm in self.s.dm:
                for ebv in self.s.ebv:
                    yield dm, ebv, make_scmd(rs, mass, apparent(absmag, dm, ebv),
                                             self.quota, self.red, self.blue)
        return age_yr, cells()

    def stage2(self, cat: Optional[np.ndarray]):
        """(JS distance, n) on the fixed grid, or None where the original skipped."""
        if cat is None:
            return None
        pts = np.column_stack([cat[:, 1] - cat[:, 2], cat[:, 0]])
        pts = pts[_window(pts, self.s.kde_window)]
        if len(pts) < 100:
            return None
        js = js_fixed(self.kde, self.p_obs, pts) if self.p_obs is not None else np.inf
        return js, len(pts)

    def stage3(self, cat: np.ndarray) -> dict:
        sim = np.column_stack([cat[:, 1] - cat[:, 2], cat[:, 0]])
        sim = sim[_window(sim, self.s.ks_window)]
        ks, d = ks_2d(self.obs3, sim, self.obs3_cdf)
        return {"ks_statistic": ks, "js_divergence": js_adaptive(self.obs3, sim),
                "mean_diff": np.mean(d), "p95_diff": np.percentile(d, 95),
                "n_observed": len(self.obs3), "n_simulated": len(sim)}


def label(age_yr: float, dm: float, ebv: float) -> Tuple[str, float, float, float]:
    """Legacy file name and the (age, DM, E(B-V)) parsed back from it."""
    name = f"NGC6791_age{age_yr / 1e9:.1f}Gyr_DM{dm:.3f}_EBV{ebv:.3f}.csv"
    return name, float(f"{age_yr / 1e9:.1f}"), float(f"{dm:.3f}"), float(f"{ebv:.3f}")


def process_age(ctx: Context, run: int, age: float, catalogs: Optional[dict] = None):
    """One (run, age): all grid cells through stage 2, the best ``keep`` through stage 3."""
    age_yr, cells = ctx.grid(run, age)
    js_rows, scored, written = [], [], []
    for dm, ebv, cat in cells:
        if cat is None:                       # stage 1 failed: no file was written
            continue
        name, a, d, e = label(age_yr, dm, ebv)
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
    """The catalogue of one grid cell, replayed from its (run, age) seed."""
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


def worker(s: Setup, out: str, runs: Sequence[int], chunk: int = 9, max_hours: float = 8.0,
           save_catalogs: bool = False, log=print) -> int:
    """Claim and process units until none are left or ``max_hours`` have passed.

    A unit is claimed by creating ``claims/<key>`` (``mkdir`` is atomic, also on
    NFS), so any number of workers -- on any account -- can share one run list
    without duplicating work.  A unit is done when its ``.ks.csv`` exists.
    """
    os.makedirs(os.path.join(out, "claims"), exist_ok=True)
    t0, n, ctx = time.time(), 0, None
    for run, ages in units(runs, s.ages, chunk):
        if time.time() - t0 > max_hours * 3600:
            break
        key = unit_key(run, ages)
        if os.path.exists(os.path.join(out, "units", str(run), key + ".ks.csv")):
            continue
        claim = os.path.join(out, "claims", key)
        try:
            os.mkdir(claim)
        except FileExistsError:
            continue
        with open(os.path.join(claim, "owner"), "w") as fh:
            fh.write(f"{os.environ.get('SLURM_JOB_ID', '-')}_{os.environ.get('SLURM_ARRAY_TASK_ID', '-')} "
                     f"{socket.gethostname()} {time.strftime('%F %T')}\n")
        if not os.path.exists(s.iso.format(run=run)):
            log(f"{key}: no isochrone file")
            continue
        ctx = ctx or Context(s)
        t = time.time()
        try:
            run_unit(ctx, out, run, ages, save_catalogs)
        except Exception:
            log(f"{key}: FAILED\n{traceback.format_exc()}")
            continue
        n += 1
        log(f"{key}: {time.time() - t:.0f} s")
    return n


def status(s: Setup, out: str, runs: Sequence[int], chunk: int = 9, release: bool = False) -> dict:
    """Units done, claimed without output yet, and untouched.

    ``release`` deletes the claims of units without output, so a new worker
    picks them up -- only do that when no worker is running.
    """
    done, claimed, todo = [], [], []
    for run, ages in units(runs, s.ages, chunk):
        key = unit_key(run, ages)
        claim = os.path.join(out, "claims", key)
        if os.path.exists(os.path.join(out, "units", str(run), key + ".ks.csv")):
            done.append(key)
        elif os.path.isdir(claim):
            owner = os.path.join(claim, "owner")
            claimed.append(f"{key} ({open(owner).read().strip() if os.path.exists(owner) else '?'})")
            if release:
                shutil.rmtree(claim)
        else:
            todo.append(key)
    return {"done": done, "claimed": claimed, "todo": todo}


def collect(s: Setup, out: str) -> pd.DataFrame:
    """Per-run ``2d_ks_results_<run>.csv`` / ``best_fit_<run>.csv`` (the legacy
    layout, so the old merge and analysis scripts read them unchanged), plus
    ``combined_2d_ks_results.csv``, ``filtered_ks_scmd.csv`` (KS <= 3 sigma)
    and every stage-2 score in ``kde_js_scores.csv.gz``."""
    ks_dir = os.path.join(out, "2d_ks")
    os.makedirs(ks_dir, exist_ok=True)
    allks, alljs = [], []
    for d in sorted(glob.glob(os.path.join(out, "units", "*"))):
        run = int(os.path.basename(d))
        ks = [pd.read_csv(f) for f in sorted(glob.glob(os.path.join(d, "*.ks.csv")))]
        js = [pd.read_csv(f) for f in sorted(glob.glob(os.path.join(d, "*.js.csv")))]
        if not ks:
            continue
        ks = pd.concat(ks, ignore_index=True)
        alljs += js
        allks.append(ks)
        ks.to_csv(os.path.join(ks_dir, f"2d_ks_results_{run}.csv"), index=False)
        best = pd.DataFrame([ks.loc[ks["ks_statistic"].idxmin()], ks.loc[ks["js_divergence"].idxmin()]])
        best["metric"] = ["best_ks", "best_js"]
        best.to_csv(os.path.join(ks_dir, f"best_fit_{run}.csv"), index=False)
    if not allks:
        return pd.DataFrame(columns=KS_COLUMNS)
    comb = pd.concat(allks, ignore_index=True).sort_values(["simulation_id", "ks_statistic"])
    comb.to_csv(os.path.join(out, "combined_2d_ks_results.csv"), index=False)
    cols = ["ks_statistic", "simulation_id", "age_gyr", "distance_modulus", "extinction_ebv"]
    comb.loc[comb["ks_statistic"] <= s.ks_3sigma, cols].to_csv(
        os.path.join(out, "filtered_ks_scmd.csv"), index=False)
    if alljs:
        pd.concat(alljs, ignore_index=True).to_csv(os.path.join(out, "kde_js_scores.csv.gz"), index=False)
    return comb
