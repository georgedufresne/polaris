"""Run the old three-script pipeline and this package side by side and compare.

The old scripts are imported unmodified; only the global NumPy RNG is seeded
(with the seed this package uses for the same (run, age)) and the per-file
loop of stage 3 is called directly so nothing is written outside a temporary
directory.  Checks, for every grid cell of the chosen (run, age):

  stage 1  the catalogue (G, BP, RP, BP-RP) -- must be bit-identical
  stage 2  the KDE Jensen-Shannon distance and the four sCMDs kept (or, when no
           sCMD qualifies, that every file of the age goes on to stage 3)
  stage 3  KS statistic, JS divergence, mean difference, star counts

Stages 2 and 3 of the old pipeline read the stage-1 CSVs back with pandas'
default float parser, which can differ from the in-memory value by one ulp,
so those are compared to a relative 1e-10 (they agree far better).

  python tests/scmd/compare_with_legacy.py --run 20001 --ages 8.3 10.5 \
      [--legacy-dir ~] [--work /tmp/x]
"""

from __future__ import annotations

import argparse
import contextlib
import glob
import importlib.util
import io
import os
import sys
import tempfile
import time

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)

from isochronetoolbox.scmd import Context, Setup, process_age, seed  # noqa: E402

STAGE1 = "scmd_grid_red_dist_chaboyer_g_err_3err_test.py"
STAGE2 = "kde_code_for_grid_best_match_chaboyer_test.py"
STAGE3 = "2d_ks_comparison_grid_chaboyer_kde.py"


def load(path: str, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def rel(a, b) -> float:
    a, b = np.asarray(a, float), np.asarray(b, float)
    fin = np.isfinite(a) & np.isfinite(b)
    if not (np.isfinite(a) == np.isfinite(b)).all():
        return np.inf
    return float(np.max(np.abs(a[fin] - b[fin]) / np.maximum(np.abs(b[fin]), 1e-300), initial=0.0))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", type=int, default=20001)
    ap.add_argument("--ages", type=float, nargs="+", default=[8.3, 10.5])
    ap.add_argument("--legacy-dir", default=os.path.expanduser("~"))
    ap.add_argument("--work", default=None)
    a = ap.parse_args()

    s = Setup()
    o1 = load(os.path.join(a.legacy_dir, STAGE1), "legacy_stage1")
    o2 = load(os.path.join(a.legacy_dir, STAGE2), "legacy_stage2")
    o3 = load(os.path.join(a.legacy_dir, STAGE3), "legacy_stage3")
    o1.DEBUG = False
    iso = s.iso.format(run=a.run)
    work = a.work or tempfile.mkdtemp(prefix="scmd_legacy_")
    grid_dir = os.path.join(work, str(a.run))
    ctx = Context(s)
    ok = True

    for age in a.ages:
        print(f"\n=== run {a.run}, age {age} Gyr ===", flush=True)
        # ---- old stage 1, seeded -------------------------------------------
        t = time.time()
        with contextlib.redirect_stdout(io.StringIO()):
            ages = o1.get_available_ages_from_isochrone(iso)
            age_yr = ages[np.argmin(np.abs(ages - age * 1e9))]
            red, blue = o1.load_regional_scatter_data(s.scatter_red, s.scatter_blue)
            np.random.seed(seed(a.run, age))
            for p in o1.create_parameter_grid_specific():
                o1.process_single_isochrone_grid_quota(
                    age_yr, p["apparent_distance_modulus"], p["ebv"], iso, s.obs, red, blue,
                    s.completeness, output_dir=grid_dir)
        t1 = time.time() - t

        # ---- new stage 1: every catalogue -----------------------------------
        _, cells = ctx.grid(a.run, age)
        n_same = n_cells = 0
        for dm, ebv, cat in cells:
            n_cells += 1
            f = os.path.join(grid_dir, o1.generate_filename(age_yr, dm, ebv, "csv"))
            old = pd.read_csv(f, float_precision="round_trip").to_numpy() if os.path.exists(f) else None
            new = None if cat is None else np.column_stack([cat, cat[:, 1] - cat[:, 2]])
            same = (old is None and new is None) or (
                old is not None and new is not None and old.shape == new.shape and np.array_equal(old, new))
            n_same += same
        print(f"stage 1: {n_same}/{n_cells} catalogues bit-identical   (old {t1:.0f} s)", flush=True)
        ok &= n_same == n_cells

        # ---- old stage 2 ------------------------------------------------------
        t = time.time()
        with contextlib.redirect_stdout(io.StringIO()):
            an = o2.GridCMDAnalyzer(str(a.run), target_age=age, grid_path=grid_dir,
                                    output_path=os.path.join(work, "kde"), obs_path=s.obs)
            res = an.run_analysis()
        t2 = time.time() - t

        # ---- new stages 2 + 3 for this age ------------------------------------
        t = time.time()
        js_new, ks_new = process_age(ctx, a.run, age)
        tn = time.time() - t
        if len(res):
            jn = pd.DataFrame(js_new).set_index(["distance_modulus", "extinction_ebv"])
            jo = res.set_index(["dm", "ebv"]).loc[jn.index]
            e2 = rel(jn["js_distance"].values, jo["js_distance"].values)
            same_n = np.array_equal(jn["n_stars"].values, jo["n_stars"].values)
            top_old = set(res.head(s.keep)["filename"])
            top_new = {r["filename"] for r in ks_new}
            print(f"stage 2: {len(jn)}/{len(res)} cells scored, JS max rel diff {e2:.1e}, "
                  f"n_stars identical {same_n}, top-{s.keep} identical {top_old == top_new}   "
                  f"(old {t2:.0f} s)", flush=True)
            ok &= len(jn) == len(res) and e2 < 1e-10 and same_n and top_old == top_new
            with contextlib.redirect_stdout(io.StringIO()):
                an.delete_non_best_files(res, top_n_to_keep=s.keep)
        else:
            # no sCMD had 100 stars in the window: the old stage 2 deleted nothing
            print(f"stage 2: no cell qualified (old and new: {len(res)} / {len(js_new)} scored); "
                  f"every file goes to stage 3   (old {t2:.0f} s)", flush=True)
            ok &= len(js_new) == 0

        # ---- old stage 3 on the survivors ------------------------------------------
        t = time.time()
        rows = []
        for f in sorted(glob.glob(os.path.join(grid_dir, f"NGC6791_age{age:.1f}Gyr_*.csv"))):
            try:                                   # as process_single_simulation does
                with contextlib.redirect_stdout(io.StringIO()):
                    r = o3.CMD_2DKS_Comparison(s.obs, f).run_comparison(
                        color_min=0.9, color_max=1.7, mag_min=16.25, mag_max=18.0)
            except Exception:
                continue
            r["filename"] = os.path.basename(f)
            rows.append(r)
        t3 = time.time() - t
        ko = pd.DataFrame(rows).set_index("filename").sort_index()
        kn = pd.DataFrame(ks_new).set_index("filename").sort_index()
        same_files = list(ko.index) == list(kn.index)
        print(f"stage 3: {len(kn)} rows new, {len(ko)} old, same files {same_files}", flush=True)
        ok &= same_files
        if not same_files:
            continue
        errs = {c: rel(kn[c].values, ko[c].values) for c in
                ("ks_statistic", "js_divergence", "mean_diff", "p95_diff")}
        cnt = all(np.array_equal(kn[c].values, ko[c].values) for c in ("n_observed", "n_simulated"))
        print("stage 3: " + ", ".join(f"{c} {v:.1e}" for c, v in errs.items())
              + f", counts identical {cnt}   (old {t3:.0f} s)", flush=True)
        ok &= all(v < 1e-10 for v in errs.values()) and cnt
        print(f"time per age: old {t1 + t2 + t3:.0f} s (stage-3 plots excluded), new {tn:.0f} s", flush=True)

    print("\nRESULT:", "IDENTICAL" if ok else "DIFFERENT", f"(work dir {work})")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
