"""Age, [Fe/H], distance and reddening from a finished grid, weighted by the bootstrap null.

The sCMD-only part of G. Dufresne's ``analyze_scmd_weighted_posteriors.py``
(NGC 6791, March-May 2026), unchanged:

1. keep the grid's stage-3 rows with KS <= mean + 3 sigma of the bootstrap
   null (``bootstrap.py``; the ``filtered_ks_scmd.csv`` that ``collect`` writes
   with the preset's ``ks_3sigma``);
2. give each row the density of the bootstrap KS distribution at its KS
   (``scipy.stats.gaussian_kde``, Scott's bandwidth);
3. keep, per (Monte Carlo run, age), the (distance, reddening) cell with the
   lowest KS;
4. normalise those densities to weights;
5. weighted histograms of age, distance modulus and reddening (from the grid)
   and of the run's Monte Carlo inputs ([Fe/H], Y, [alpha/Fe], mixing length),
   each summarised by mean, median, mode, 16th/84th percentiles (interpolated
   in the cumulative histogram) and standard deviation (``compute_stats``).

With NGC 6791's inputs (``filtered_ks_scmd.csv``, ``Bootstrap/
bootstrap_2d_ks_results.csv``, ``initial_vars.csv``) this gives the original's
numbers (``tests/scmd/test_posterior.py``).
"""

from __future__ import annotations

import glob
import os
import re
from typing import Dict, Optional, Sequence

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde


def _edges(start, stop, step):
    e = np.arange(start, stop, step)
    return e, (e[:-1] + e[1:]) / 2


def _bins_ngc6791() -> Dict[str, dict]:
    """The original's bins."""
    b = {"age": _edges(3.95, 12.15, 0.1), "feh": _edges(-0.505, 1.005, 0.02),
         "dm": _edges(10.005, 16.005, 0.02), "ebv": _edges(-0.005, 0.805, 0.01),
         "Y": _edges(0.195, 0.405, 0.002), "cmix": _edges(0.995, 2.505, 0.01)}
    out = {k: {"edges": e, "centers": c} for k, (e, c) in b.items()}
    out["afe"] = {"edges": np.array([-0.3, -0.1, 0.1, 0.3, 0.5]), "centers": np.array([-0.2, 0.0, 0.2, 0.4]),
                  "discrete": True}
    return out


def _bins_47tuc() -> Dict[str, dict]:
    """Bins on the 47 Tuc grid's steps: 0.1 Gyr, mu 0.01, E(V-I) 0.005."""
    b = {"age": _edges(7.95, 14.15, 0.1), "feh": _edges(-1.205, -0.295, 0.02),
         "dm": _edges(13.145, 13.465, 0.01), "ebv": _edges(-0.0025, 0.0551, 0.005),
         "Y": _edges(0.195, 0.405, 0.002), "cmix": _edges(0.995, 2.505, 0.01)}
    out = {k: {"edges": e, "centers": c} for k, (e, c) in b.items()}
    out["afe"] = {"edges": np.array([-0.1, 0.1, 0.3, 0.5]), "centers": np.array([0.0, 0.2, 0.4]),
                  "discrete": True}
    return out


BINS = {"ngc6791": _bins_ngc6791, "47tuc": _bins_47tuc}

# Monte Carlo inputs per setup name
MC_INPUTS = {"NGC6791": os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
                                     "data", "ngc6791", "mc_initial_vars.csv"),
             "47Tuc": "/dartfs-hpc/scratch/f007058/47Tuc_var/varfeh-077.{run}"}

# parameter -> (column in the merged table, label)
PARAMS = {"age": ("age_gyr", "age (Gyr)"), "feh": ("FeH", "[Fe/H]"),
          "dm": ("distance_modulus", "distance modulus"), "ebv": ("extinction_ebv", "reddening"),
          "Y": ("Y", "Y"), "afe": ("alphafe", "[alpha/Fe]"), "cmix": ("CMIXLA", "mixing length")}


def compute_stats(vals, prob, discrete=False) -> dict:
    """Mean, median, mode, 16th/84th percentiles and std of a histogram (the original's)."""
    prob = prob / np.sum(prob)
    mean = np.sum(vals * prob)
    cum = np.cumsum(prob)

    def interp_pct(cum, vals, pct):
        idx = np.searchsorted(cum, pct)
        if idx == 0:
            return vals[0]
        if idx >= len(vals):
            return vals[-1]
        if discrete:
            return vals[idx]
        frac = (pct - cum[idx - 1]) / (cum[idx] - cum[idx - 1] + 1e-10)
        return vals[idx - 1] + frac * (vals[idx] - vals[idx - 1])

    return {"mean": mean, "median": interp_pct(cum, vals, 0.50), "mode": vals[np.argmax(prob)],
            "p16": interp_pct(cum, vals, 0.16), "p84": interp_pct(cum, vals, 0.84),
            "std": np.sqrt(np.sum(prob * (vals - mean) ** 2))}


def read_mc_inputs(source: str, runs: Optional[Sequence[int]] = None) -> pd.DataFrame:
    """Monte Carlo inputs per run, with ``simulation_id`` and ``Y`` columns.

    ``source`` is a CSV with ``file_id`` (NGC 6791's ``initial_vars.csv``) or a
    pattern with ``{run}`` for DSEP ``varfeh`` files ("value name" per line;
    47 Tuc), read for ``runs`` (default: every file matching the pattern)."""
    if "{run}" not in source:
        df = pd.read_csv(source)
        df["simulation_id"] = df.pop("file_id").astype(int)
    else:
        if runs is None:
            rx = re.compile(re.escape(source).replace(r"\{run\}", r"(\d+)") + "$")
            runs = sorted(int(rx.match(f).group(1)) for f in glob.glob(source.replace("{run}", "*"))
                          if rx.match(f))
        rec = []
        for r in runs:
            path = source.format(run=r)
            if not os.path.exists(path):
                continue
            v = {}
            for line in open(path):
                parts = line.split(None, 1)
                if len(parts) == 2:
                    v[parts[1].strip()] = float(parts[0])
            v["simulation_id"] = int(r)
            rec.append(v)
        df = pd.DataFrame(rec)
    df = df.rename(columns={"He Abundance": "Y", "Yprim": "Y"})
    return df


def weigh(rows: pd.DataFrame, boot_ks: np.ndarray, limit: Optional[float] = None,
          cut: bool = True) -> pd.DataFrame:
    """Steps 1-4: the rows that carry weight, with ``KS_Probability_Density`` and
    ``P_KS_normalized``; empty if no row passes the cut.  ``limit`` defaults to
    the bootstrap's mean + 3 sigma (NGC 6791's grid was cut at it rounded, 0.0468)."""
    boot_ks = np.asarray(boot_ks, float)
    if limit is None:
        limit = np.mean(boot_ks) + 3 * np.std(boot_ks)
    df = rows[rows["ks_statistic"] <= limit].copy() if cut else rows.copy()
    if df.empty:
        return df.assign(KS_Probability_Density=[], P_KS_normalized=[])
    df["KS_Probability_Density"] = gaussian_kde(boot_ks, bw_method="scott")(df["ks_statistic"].values)
    df = df.loc[df.groupby(["simulation_id", "age_gyr"])["ks_statistic"].idxmin()]
    total = df["KS_Probability_Density"].sum()
    df["P_KS_normalized"] = df["KS_Probability_Density"] / total if total > 0 else np.nan
    return df


def posteriors(weighted: pd.DataFrame, mc: Optional[pd.DataFrame], bins: Dict[str, dict]):
    """Step 5: {parameter: (centers, probabilities, stats)} for every parameter with data."""
    df = weighted.merge(mc, on="simulation_id", how="left") if mc is not None else weighted
    out = {}
    for key, (col, _) in PARAMS.items():
        if col not in df.columns or key not in bins:
            continue
        vals = df[col].dropna()
        if vals.empty:
            continue
        w = df.loc[vals.index, "P_KS_normalized"].values
        b = bins[key]
        h, _ = np.histogram(vals.values, bins=b["edges"], weights=w)
        if h.sum() > 0:
            h = h / h.sum()
            out[key] = (b["centers"], h, compute_stats(b["centers"], h, b.get("discrete", False)))
    return out


def summary_table(post: dict, weighted: pd.DataFrame) -> pd.DataFrame:
    rec = []
    for key, (_, _, st) in post.items():
        rec.append({"parameter": key, **st, "rows": len(weighted),
                    "runs": int(weighted["simulation_id"].nunique())})
    return pd.DataFrame(rec)


def plot(post: dict, path: str, title: str, reddening_label: str = "E(B-V)") -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    keys = [k for k in PARAMS if k in post]
    fig, axes = plt.subplots(1, len(keys), figsize=(3.6 * len(keys), 3.4))
    for ax, k in zip(np.atleast_1d(axes), keys):
        c, h, st = post[k]
        width = np.median(np.diff(c)) if len(c) > 1 else 0.1
        ax.bar(c, h, width=width * 0.92, color="#4a6fa5")
        ax.axvline(st["median"], color="0.15", lw=1.2)
        ax.axvspan(st["p16"], st["p84"], color="0.6", alpha=0.25, lw=0)
        lab = reddening_label if k == "ebv" else PARAMS[k][1]
        ax.set_xlabel(lab)
        ax.set_title(f"{st['median']:.3f} (+{st['p84'] - st['median']:.3f} / -{st['median'] - st['p16']:.3f})",
                     fontsize=9)
        nz = c[h > 0]
        if len(nz):
            pad = 3 * width
            ax.set_xlim(nz.min() - pad, nz.max() + pad)
    np.atleast_1d(axes)[0].set_ylabel("posterior probability")
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)
