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

from .pipeline import Context, Setup, collect, regenerate, run_unit, status, units, worker


def _setup(a) -> Setup:
    path = os.path.join(a.out, "setup.json")
    given = {k: v for k, v in (("iso", a.iso), ("obs", a.obs), ("bc", a.bc)) if v}
    if a.exact_age:
        given["exact_age"] = True
    if os.path.exists(path):
        saved = json.load(open(path))
        clash = {k: (v, saved.get(k)) for k, v in given.items() if saved.get(k) != v}
        if clash:
            sys.exit(f"{path} was made with a different setup: {clash}")
        given = saved
    s = Setup(**{k: v for k, v in given.items() if k in ("iso", "obs", "bc", "exact_age")})
    if not os.path.exists(path):
        os.makedirs(a.out, exist_ok=True)
        keep = {k: v for k, v in dataclasses.asdict(s).items() if isinstance(v, (str, bool, int, float))}
        tmp = f"{path}.tmp{os.getpid()}"
        json.dump(keep, open(tmp, "w"), indent=1)
        os.replace(tmp, path) if not os.path.exists(path) else os.remove(tmp)
    return s


def _runs(path: str):
    return [int(x) for x in open(path).read().split()]


def main(argv=None) -> None:
    p = argparse.ArgumentParser(prog="python -m isochronetoolbox.scmd")
    p.add_argument("command", choices=["worker", "run", "status", "collect", "plot"])
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--runs", help="file of run numbers (worker, status)")
    p.add_argument("--run", type=int, help="one run number (run)")
    p.add_argument("--chunk", type=int, default=9, help="ages per unit of work")
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
        n = worker(s, a.out, _runs(a.runs), a.chunk, a.max_hours, a.save_catalogs,
                   log=lambda m: print(m, flush=True))
        print(f"worker done: {n} units", flush=True)
    elif a.command == "run":
        ctx = Context(s)
        for run, ages in units([a.run], s.ages, a.chunk):
            run_unit(ctx, a.out, run, ages, a.save_catalogs)
    elif a.command == "status":
        st = status(s, a.out, _runs(a.runs), a.chunk, a.release)
        print({k: len(v) for k, v in st.items()})
        if st["claimed"]:
            print(("released:" if a.release else "claimed, no output yet:"), *st["claimed"][:40], sep="\n  ")
    elif a.command == "collect":
        comb = collect(s, a.out)
        print(f"{comb['simulation_id'].nunique()} runs, {len(comb)} rows, "
              f"{(comb['ks_statistic'] <= s.ks_3sigma).sum()} with KS <= {s.ks_3sigma}")
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
            cat = regenerate(ctx, r.simulation_id, r.age_gyr, r.distance_modulus, r.extinction_ebv)
            title = (f"#{rank}  run {r.simulation_id}   age {r.age_gyr:.1f} Gyr   DM {r.distance_modulus:.2f}"
                     f"   E(B-V) {r.extinction_ebv:.2f}   KS {r.ks_statistic:.4f}   JS {r.js_divergence:.4f}")
            path = os.path.join(d, f"rank{rank:02d}_run{r.simulation_id}_{r.filename.replace('.csv', '.png')}")
            cmd_pair(obs, cat, title, path, ks_window=s.ks_window)
            print(path)


if __name__ == "__main__":
    main()
