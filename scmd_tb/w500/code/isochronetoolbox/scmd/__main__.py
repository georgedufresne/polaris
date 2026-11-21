"""Command line: ``python -m isochronetoolbox.scmd <command> --out DIR ...``

  worker   claim units from a run list and process them (one per SLURM task)
  run      process one run directly, no claims
  status   done / claimed / to-do units for a run list
  collect  legacy per-run 2D-KS tables, combined table, 3-sigma filtered table
  plot     observed CMD next to the best-fitting sCMDs (regenerated from seed)

The grid definition is written to ``DIR/setup.json`` by the first command
that uses ``DIR`` and read back by every later one, so one output directory
can never mix two configurations.
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import sys

import numpy as np
import pandas as pd

from .pipeline import Context, Setup, collect, regenerate, run_unit, status, unit_key, units, worker


_SAVED = ("iso", "obs", "completeness", "scatter_red", "scatter_blue", "bc",
          "exact_age", "keep", "chunk", "ks_3sigma")


def _setup(a) -> Setup:
    """The grid definition of ``a.out``: created atomically by the first command,
    and every later one -- or any concurrent worker -- must agree with it."""
    path = os.path.join(a.out, "setup.json")
    given = {k: os.path.abspath(getattr(a, k)) for k in ("iso", "obs", "bc") if getattr(a, k)}
    if a.exact_age:
        given["exact_age"] = True
    if a.chunk is not None:
        given["chunk"] = a.chunk
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
    clash = {k: (v, saved.get(k, getattr(Setup, k, None))) for k, v in given.items()
             if saved.get(k, getattr(Setup, k, None)) != v}
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
    p.add_argument("command", choices=["worker", "run", "status", "collect", "plot"])
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--runs", help="file of run numbers (worker, status)")
    p.add_argument("--run", type=int, help="one run number (run)")
    p.add_argument("--chunk", type=int, help="ages per unit of work (default 9, fixed per --out)")
    p.add_argument("--max-hours", type=float, default=8.0, help="stop claiming after this")
    p.add_argument("--save-catalogs", action="store_true", help="also keep the stage-3 sCMDs")
    p.add_argument("--iso", help="isochrone path pattern with {run}")
    p.add_argument("--obs", help="observed catalogue (BP_RP or BPmag/RPmag, Gmag)")
    p.add_argument("--bc", help="BC table (.npz); default MIST [Fe/H] = -3.0 as before")
    p.add_argument("--exact-age", action="store_true", help="no 1%% age union above 10 Gyr")
    p.add_argument("--top", type=int, default=10, help="plot: how many sCMDs")
    p.add_argument("--distinct-runs", action="store_true", help="plot: at most one per run")
    p.add_argument("--release", action="store_true",
                   help="status: free the claims of units without output (no workers may be running)")
    a = p.parse_args(argv)
    s = _setup(a)

    if a.command == "worker":
        n = worker(s, a.out, _runs(a.runs), a.max_hours, a.save_catalogs,
                   log=lambda m: print(m, flush=True))
        print(f"worker done: {n} units", flush=True)
    elif a.command == "run":
        ctx = Context(s)
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
    elif a.command == "plot":
        from .plot import cmd_pair
        comb = pd.read_csv(os.path.join(a.out, "combined_2d_ks_results.csv")).sort_values("ks_statistic")
        if a.distinct_runs:
            comb = comb.drop_duplicates("simulation_id")
        best = comb.head(a.top)
        ctx, obs = Context(s), pd.read_csv(s.obs)
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
            cmd_pair(obs, cat, title, path, ks_window=s.ks_window)
            print(path)


if __name__ == "__main__":
    main()
