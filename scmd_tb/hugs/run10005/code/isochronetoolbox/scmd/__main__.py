"""Command line: ``python -m isochronetoolbox.scmd <command> --out DIR ...``

  worker   claim units from a run list and process them (one per SLURM task)
  run      process one run directly, no claims
  status   done / claimed / to-do units for a run list
  collect  legacy per-run 2D-KS tables, combined table, 3-sigma filtered table
  plot     observed CMD next to the best-fitting sCMDs (regenerated from seed)
  bootstrap          realisations of one cell scored against its grid sCMD (bootstrap.py);
                     a slice per SLURM task, parts in DIR/bootstrap/
  bootstrap-collect  merge the parts: DIR/bootstrap_2d_ks_results.csv, bootstrap_summary.json
  posterior          weighted age, [Fe/H], distance, reddening, ... (posterior.py) in DIR/posterior/

The grid definition is written to ``DIR/setup.json`` by the first command
that uses ``DIR`` and read back by every later one, so one output directory
can never mix two configurations.
"""

from __future__ import annotations

import argparse
import dataclasses
import filecmp
import json
import os
import sys

import numpy as np
import pandas as pd

from .pipeline import (PRESETS, Setup, collect, context, regenerate, run_unit, status, unit_key,
                       units, worker)


_SAVED = tuple(f.name for f in dataclasses.fields(Setup))


def _norm(v):
    return json.loads(json.dumps(v))


def _same(a, b) -> bool:
    """Setup values agree; two paths to identical files (e.g. the package data
    seen from a frozen copy of the code) count as the same."""
    if a == b:
        return True
    if isinstance(a, str) and isinstance(b, str) and os.path.isfile(a) and os.path.isfile(b):
        return filecmp.cmp(a, b, shallow=False)
    return False


def _setup(a) -> Setup:
    """The grid definition of ``a.out``: created atomically by the first command,
    and every later one -- or any concurrent worker -- must agree with it."""
    path = os.path.join(a.out, "setup.json")
    given = dict(PRESETS[a.preset.lower()]) if a.preset else {}
    given.update({k: os.path.abspath(getattr(a, k)) for k in ("iso", "obs", "bc") if getattr(a, k)})
    if a.exact_age:
        given["exact_age"] = True
    if a.chunk is not None:
        given["chunk"] = a.chunk
    if a.ridge_width is not None:
        given["ridge_width"] = a.ridge_width
    if a.ks_percentile is not None:
        given["ks_percentile"] = a.ks_percentile
    os.makedirs(a.out, exist_ok=True)
    if not os.path.exists(path):
        rec = {k: v for k, v in dataclasses.asdict(Setup(**given)).items() if k in _SAVED}
        tmp = f"{path}.tmp{os.getpid()}"
        with open(tmp, "w") as fh:
            json.dump(rec, fh, indent=1)
        try:
            os.link(tmp, path)            # atomic: exactly one process creates it
        except FileExistsError:
            pass
        finally:
            os.remove(tmp)
    saved = json.load(open(path))
    default = {k: _norm(v) for k, v in dataclasses.asdict(Setup()).items()}
    clash = {k: (_norm(v), saved.get(k, default.get(k))) for k, v in given.items()
             if not _same(_norm(v), saved.get(k, default.get(k)))}
    if clash:
        sys.exit(f"{path} was made with a different setup: {clash}")
    return Setup(**{k: v for k, v in saved.items() if k in _SAVED})


def _saved(out: str, s: Setup, r):
    """The scored catalogue of a table row, if the worker kept it (--save-catalogs)."""
    for run, ages in units([r.simulation_id], s.ages, s.chunk):
        if any(abs(age - r.age_gyr) < 1e-6 for age in ages):
            path = os.path.join(out, "units", str(run), unit_key(run, ages) + ".npz")
            if os.path.exists(path):
                z = np.load(path)
                if r.filename in z.files:
                    return z[r.filename]
    return None


def _runs(path: str):
    return [int(x) for x in open(path).read().split()]


def main(argv=None) -> None:
    p = argparse.ArgumentParser(prog="python -m isochronetoolbox.scmd")
    p.add_argument("command", choices=["worker", "run", "status", "collect", "plot", "bootstrap",
                                       "bootstrap-collect", "posterior"])
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--preset", help="named setup: " + ", ".join(PRESETS) + " (default ngc6791)")
    p.add_argument("--runs", help="file of run numbers (worker, status)")
    p.add_argument("--run", type=int, help="one run number (run)")
    p.add_argument("--chunk", type=int, help="ages per unit of work (default 9, fixed per --out)")
    p.add_argument("--max-hours", type=float, default=8.0, help="stop claiming after this")
    p.add_argument("--save-catalogs", action="store_true", help="also keep the stage-3 sCMDs")
    p.add_argument("--iso", help="isochrone path pattern with {run}")
    p.add_argument("--obs", help="observed catalogue (BP_RP or BPmag/RPmag, Gmag)")
    p.add_argument("--bc", help="BC table (.npz); default MIST [Fe/H] = -3.0 as before")
    p.add_argument("--exact-age", action="store_true", help="no 1%% age union above 10 Gyr")
    p.add_argument("--ridge-width", type=float, help="cut=ridge: normalised perpendicular distance kept")
    p.add_argument("--ks-percentile", type=float, help="stage-3 2D KS percentile of |ECDF diff| (default 95; fixed per --out)")
    p.add_argument("--top", type=int, default=10, help="plot: how many sCMDs")
    p.add_argument("--distinct-runs", action="store_true", help="plot: at most one per run")
    p.add_argument("--release", action="store_true",
                   help="status: free the claims of units without output (no workers may be running)")
    p.add_argument("--cell", help="bootstrap: RUN:AGE:DM:EBV (default: the lowest KS in combined_2d_ks_results.csv)")
    p.add_argument("--n", type=int, default=10000, help="bootstrap: number of realisations (0 .. n-1)")
    p.add_argument("--ntasks", type=int, help="bootstrap: total number of slices over all arrays "
                   "(default SLURM_ARRAY_TASK_COUNT, else 1)")
    p.add_argument("--task-offset", type=int, default=0,
                   help="bootstrap: added to SLURM_ARRAY_TASK_ID, so arrays on two accounts take distinct slices")
    p.add_argument("--task", type=int, help="bootstrap: this slice (default offset + SLURM_ARRAY_TASK_ID, else 0)")
    p.add_argument("--bootstrap", help="posterior: bootstrap KS table (default DIR/bootstrap_2d_ks_results.csv)")
    p.add_argument("--mc", help="posterior: Monte Carlo inputs, a CSV with file_id or a pattern with {run}")
    p.add_argument("--limit", type=float, help="posterior: KS cut (default bootstrap mean + 3 sigma)")
    a = p.parse_args(argv)
    s = _setup(a)

    if a.command == "worker":
        n = worker(s, a.out, _runs(a.runs), a.max_hours, a.save_catalogs,
                   log=lambda m: print(m, flush=True))
        print(f"worker done: {n} units", flush=True)
    elif a.command == "run":
        ctx = context(s)
        for run, ages in units([a.run], s.ages, s.chunk):
            run_unit(ctx, a.out, run, ages, a.save_catalogs)
    elif a.command == "status":
        st = status(s, a.out, _runs(a.runs), a.release)
        print({k: len(v) for k, v in st.items()})
        if st["claimed"]:
            print(("released:" if a.release else "claimed, no output yet:"), *st["claimed"][:40], sep="\n  ")
    elif a.command == "collect":
        comb, rep = collect(s, a.out)
        print(f"{len(rep['complete'])} complete runs ({len(rep['no_rows'])} without stage-3 rows), "
              f"{len(comb)} rows, {(comb['ks_statistic'] <= s.ks_3sigma).sum()} with KS <= {s.ks_3sigma}")
        if rep["incomplete"]:
            print(f"left out, units missing: {len(rep['incomplete'])} runs:", *rep["incomplete"][:30])
    elif a.command == "bootstrap":
        from . import bootstrap as B
        if a.cell:
            r, ag, dm, e = a.cell.split(":")
            cell = {"simulation_id": int(r), "age_gyr": float(ag), "distance_modulus": float(dm),
                    "extinction_ebv": float(e)}
        else:
            comb = pd.read_csv(os.path.join(a.out, "combined_2d_ks_results.csv"))
            cell = comb.loc[comb["ks_statistic"].idxmin()].to_dict()
        env = os.environ
        ntasks = a.ntasks or int(env.get("SLURM_ARRAY_TASK_COUNT", 1))
        task = a.task if a.task is not None else a.task_offset + int(env.get("SLURM_ARRAY_TASK_ID", 0))
        part = B.run_slice(context(s), a.out, cell, a.n, task, ntasks, log=lambda m: print(m, flush=True))
        print(f"bootstrap slice {task}/{ntasks} -> {part}", flush=True)
    elif a.command == "bootstrap-collect":
        from . import bootstrap as B
        rep = B.collect(a.out)
        print(f"{rep['n'] - rep['missing']} of {rep['n']} realisations: KS mean {rep['ks_mean']:.5f}, "
              f"std {rep['ks_std']:.5f}, mean + 3 sigma {rep['ks_3sigma']:.5f}")
    elif a.command == "posterior":
        from . import posterior as P
        rows = pd.read_csv(os.path.join(a.out, "combined_2d_ks_results.csv"))
        boot = pd.read_csv(a.bootstrap or os.path.join(a.out, "bootstrap_2d_ks_results.csv"))["ks_statistic"].values
        mc_src = a.mc or P.MC_INPUTS.get(s.name)
        mc = P.read_mc_inputs(mc_src, sorted(rows["simulation_id"].unique())) if mc_src else None
        limit = a.limit if a.limit is not None else float(np.mean(boot) + 3 * np.std(boot))
        w = P.weigh(rows, boot, limit)
        print(f"bootstrap: {len(boot)} realisations, KS mean {np.mean(boot):.5f}, std {np.std(boot):.5f}; "
              f"cut KS <= {limit:.5f}: {(rows['ks_statistic'] <= limit).sum()} of {len(rows)} rows; "
              f"{len(w)} (run, age) best cells from {w['simulation_id'].nunique() if len(w) else 0} runs carry weight")
        if w.empty:
            print(f"no sCMD has KS <= {limit:.5f} (best {rows['ks_statistic'].min():.5f}): no posterior")
            return
        post = P.posteriors(w, mc, P.BINS.get(s.name.lower(), P.BINS["ngc6791"])())
        d = os.path.join(a.out, "posterior")
        os.makedirs(d, exist_ok=True)
        tab = P.summary_table(post, w)
        tab.to_csv(os.path.join(d, "posterior_summary.csv"), index=False)
        w.to_csv(os.path.join(d, "posterior_weights.csv"), index=False)
        for k, (c, h, _) in post.items():
            pd.DataFrame({"center": c, "probability": h}).to_csv(os.path.join(d, f"hist_{k}.csv"), index=False)
        red = "E(V-I)" if s.reddening == "evi" else "E(B-V)"
        P.plot(post, os.path.join(d, "posteriors.png"), f"{s.name}: sCMD posteriors "
               f"({len(w)} best cells, {w['simulation_id'].nunique()} runs, KS <= {limit:.4f})", red)
        print(tab.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    elif a.command == "plot":
        from .plot import cmd_pair
        comb = pd.read_csv(os.path.join(a.out, "combined_2d_ks_results.csv")).sort_values("ks_statistic")
        if a.distinct_runs:
            comb = comb.drop_duplicates("simulation_id")
        best = comb.head(a.top)
        ctx = context(s)
        d = os.path.join(a.out, "plots")
        os.makedirs(d, exist_ok=True)
        for rank, r in enumerate(best.itertuples(), 1):
            cat = _saved(a.out, s, r)
            if cat is None:
                # replayed from the seed; isochrone sets with tied masses replay
                # differently on another CPU type (AVX-512 vs AVX2), so check it
                cat = regenerate(ctx, r.simulation_id, r.age_gyr, r.distance_modulus, r.extinction_ebv)
                ks = ctx.stage3(cat)["ks_statistic"]
                if abs(ks - r.ks_statistic) > 1e-9:
                    print(f"WARNING rank {rank}: the replayed sCMD scores KS {ks:.4f}, the table {r.ks_statistic:.4f};"
                          " it was scored on another CPU type -- plot on the workers' node type or use --save-catalogs")
            title = (f"#{rank}  run {r.simulation_id}   age {r.age_gyr:.1f} Gyr   DM {r.distance_modulus:.2f}"
                     f"   E(B-V) {r.extinction_ebv:.2f}   KS {r.ks_statistic:.4f}   JS {r.js_divergence:.4f}")
            path = os.path.join(d, f"rank{rank:02d}_run{r.simulation_id}_{r.filename.replace('.csv', '.png')}")
            if s.model == "ast":
                (c0, c1), (m0, m1) = s.ks_window
                cmd_pair(ctx.obs, ctx.points(cat), title, path, window=s.ks_window,
                         labels=("F606W - F814W", "F606W"), xlim=(0.35, 1.05), ylim=(m1 + 0.2, m0 - 0.2))
            else:
                cmd_pair(ctx.obs, ctx.points(cat), title, path, window=s.ks_window)
            print(path)


if __name__ == "__main__":
    main()
