"""The bootstrap null of the 2D KS statistic.

How large is the KS statistic when the model is right?  G. Dufresne's NGC 6791
answer (``scmd_bootstrap.py``, ``2d_ks_bootstrap_kde_w_.py``, July 2025): take
the best-fitting grid cell (run 20140, 8.4 Gyr, DM 13.33, E(B-V) 0.16), use its
sCMD from the grid as a stand-in for the observed CMD, draw 10,000 more sCMDs
of the same cell that differ only in their random numbers (seed = bootstrap
number 0-9999), and score each against the stand-in with stage 3 (2D KS in the
KS window).  The spread of those values is what sampling alone does: NGC 6791 has
mean 0.0268, standard deviation 0.0067, so KS <= mean + 3 sigma = 0.0468 marks
a fit consistent with the data, and the posteriors weight each fit by the
density of this distribution at its KS (``posterior.py``).

Here every realisation is made by the grid's own generator (``Context.model``,
so the same model, star count and selection as the grid) and scored by the
grid's own stage 3 (``Context.compare``).  The realisations are independent of
each other: realisation ``i`` uses ``np.random.RandomState(i)``, as the
original seeded with the bootstrap number, so any subset can be computed on any
SLURM task and the union is the same.

The null depends on the star count: ~1,100 stars in the NGC 6791 KS window,
~47,000 for 47 Tuc, where it is ~6x narrower.
"""

from __future__ import annotations

import glob
import json
import os
from typing import Iterable, List

import numpy as np
import pandas as pd

from .pipeline import Context, _window, regenerate
from .stats import rank2d

COLUMNS = ["bootstrap_number", "ks_statistic", "js_divergence", "mean_diff", "p95_diff",
           "n_observed", "n_simulated", "simulation_id", "age_gyr", "distance_modulus",
           "extinction_ebv"]


def reference(ctx: Context, run: int, age: float, dm: float, ebv: float) -> np.ndarray:
    """The stand-in for the observed CMD: the grid's sCMD of this cell, cut to the KS window."""
    pts = ctx.points(regenerate(ctx, run, age, dm, ebv))
    return pts[_window(pts, ctx.s.ks_window)]


def realisations(ctx: Context, run: int, age: float, dm: float, ebv: float,
                 ids: Iterable[int], log=None) -> List[dict]:
    """Stage-3 rows of realisations ``ids`` of one cell against its grid sCMD."""
    ref = reference(ctx, run, age, dm, ebv)
    ref_cdf = rank2d(ref) / len(ref)
    _, make = ctx.model(run, age)
    rows = []
    for i in ids:
        cat = make(np.random.RandomState(int(i)), dm, ebv)
        row = ctx.compare(ref, ref_cdf, cat)
        row.update(bootstrap_number=int(i), simulation_id=int(run), age_gyr=float(age),
                   distance_modulus=float(dm), extinction_ebv=float(ebv))
        rows.append(row)
        if log:
            log(f"bootstrap {i}: KS {row['ks_statistic']:.6f}  JS {row['js_divergence']:.6f}")
    return rows


def summary(ks: np.ndarray) -> dict:
    """Mean, standard deviation (ddof = 0, as the posterior script) and the 3-sigma limit."""
    ks = np.asarray(ks, float)
    mean, std = float(np.mean(ks)), float(np.std(ks))
    return {"n": int(len(ks)), "ks_mean": mean, "ks_std": std, "ks_median": float(np.median(ks)),
            "ks_min": float(np.min(ks)), "ks_max": float(np.max(ks)), "ks_3sigma": mean + 3 * std}


def run_slice(ctx: Context, out: str, cell: dict, n: int, task: int, ntasks: int, log=None) -> str:
    """Realisations ``task, task + ntasks, ...`` below ``n`` into ``out/bootstrap/part_TASK.csv``.

    The cell is recorded in ``out/bootstrap/cell.json`` by the first task; every
    other task must agree with it, so two bootstraps can never mix."""
    d = os.path.join(out, "bootstrap")
    os.makedirs(d, exist_ok=True)
    want = {"simulation_id": int(cell["simulation_id"]), "age_gyr": float(cell["age_gyr"]),
            "distance_modulus": round(float(cell["distance_modulus"]), 4),
            "extinction_ebv": round(float(cell["extinction_ebv"]), 4), "n": int(n)}
    path = os.path.join(d, "cell.json")
    if not os.path.exists(path):
        tmp = f"{path}.tmp{os.getpid()}"
        with open(tmp, "w") as fh:
            json.dump(want, fh, indent=1)
        try:
            os.link(tmp, path)
        except FileExistsError:
            pass
        finally:
            os.remove(tmp)
    have = json.load(open(path))
    if have != want:
        raise SystemExit(f"{path} holds another bootstrap: {have} (asked for {want})")
    ids = range(task, n, ntasks)
    rows = realisations(ctx, want["simulation_id"], want["age_gyr"], want["distance_modulus"],
                        want["extinction_ebv"], ids, log)
    part = os.path.join(d, f"part_{task:05d}.csv")
    pd.DataFrame(rows)[COLUMNS].to_csv(f"{part}.tmp", index=False)
    os.replace(f"{part}.tmp", part)
    return part


def collect(out: str) -> dict:
    """Merge the parts into ``out/bootstrap_2d_ks_results.csv`` and write the summary
    to ``out/bootstrap_summary.json``; realisations missing are reported."""
    d = os.path.join(out, "bootstrap")
    cell = json.load(open(os.path.join(d, "cell.json")))
    parts = sorted(glob.glob(os.path.join(d, "part_*.csv")))
    df = pd.concat([pd.read_csv(p) for p in parts]).sort_values("bootstrap_number")
    missing = sorted(set(range(cell["n"])) - set(df["bootstrap_number"]))
    df.to_csv(os.path.join(out, "bootstrap_2d_ks_results.csv"), index=False)
    rep = {**cell, **summary(df["ks_statistic"]), "missing": len(missing)}
    with open(os.path.join(out, "bootstrap_summary.json"), "w") as fh:
        json.dump(rep, fh, indent=1)
    return rep
