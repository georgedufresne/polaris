"""The posterior port against the original sCMD-only posterior (NGC 6791).

Runs sections 1, 1b and 1c of ``analyze_scmd_weighted_posteriors.py`` on its own
inputs and compares every statistic with ``posterior.py``.  Needs the original
script and its inputs in the author's home directory; skipped elsewhere.

Run directly (``python tests/scmd/test_posterior.py``) or under pytest.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from isochronetoolbox.scmd import posterior as P  # noqa: E402

HOME = "/dartfs-hpc/rc/home/8/f007058"
SCRIPT = f"{HOME}/analyze_scmd_weighted_posteriors.py"
INPUTS = [f"{HOME}/filtered_ks_scmd.csv", f"{HOME}/Bootstrap/bootstrap_2d_ks_results.csv",
          f"{HOME}/initial_vars.csv"]


def _original():
    from scipy.stats import gaussian_kde
    src = open(SCRIPT).read()
    cfg = src[src.index("OUTPUT_DIR ="):src.index("# HELPER FUNCTION")]
    helper = src[src.index("def compute_stats"):src.index("# 1. LOAD sCMD DATA")]
    s1 = src[src.index("filtered_scmd = pd.read_csv"):src.index("# 1a-bis.")]
    s1b = src[src.index("initial_vars = pd.read_csv"):src.index("n_feh_ok =")]
    s1c = src[src.index("scmd_only_posteriors = {}"):src.index("# 2. LOAD DEB CHI2")]
    g = {"np": np, "pd": pd, "gaussian_kde": gaussian_kde, "os": os, "print": lambda *a, **k: None}
    exec(cfg.split("# 9 tests")[0] + cfg[cfg.index("# Define age grids"):], g)
    exec(helper.rsplit("# ===", 1)[0], g)
    exec(s1 + "\n" + s1b + "\n" + s1c.rsplit("# ===", 1)[0], g)
    return g["scmd_only_stats"]


def test_reproduces_original_posterior():
    if not all(os.path.exists(p) for p in [SCRIPT] + INPUTS):
        print("original posterior or its inputs not found: skipped")
        return
    ref = _original()
    rows = pd.read_csv(INPUTS[0])
    boot = pd.read_csv(INPUTS[1])["ks_statistic"].values
    w = P.weigh(rows, boot, limit=0.0468)
    post = P.posteriors(w, P.read_mc_inputs(INPUTS[2]), P.BINS["ngc6791"]())
    for k_ref, k in {"age": "age", "feh": "feh", "dm": "dm", "ebv": "ebv", "he": "Y",
                     "afe": "afe", "cmix": "cmix"}.items():
        for stat, v in ref[k_ref].items():
            assert abs(post[k][2][stat] - v) < 1e-12, (k, stat, post[k][2][stat], v)


def test_compute_stats():
    c = np.array([1.0, 2.0, 3.0])
    st = P.compute_stats(c, np.array([0.25, 0.5, 0.25]))
    assert st["mean"] == 2.0 and st["mode"] == 2.0 and abs(st["std"] - np.sqrt(0.5)) < 1e-12


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(name, "OK")
