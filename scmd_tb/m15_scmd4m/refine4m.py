#!/usr/bin/env python3
"""KS descent with 4M-star model sCMDs, from the 32k-star refined minima of a finished M15 grid.

Why: every production sCMD has as many stars as the observed fit sample (32,085), so the model ECDF is as
noisy as the data's. Ying et al. (2025) use 4,000,000 points per isochrone. This rescoring keeps everything
else (setup, ridge cut traced on an observed-sized sample, nearest-AST errors, crn replay, 95th-percentile
2D KS) and changes only the model size (``n_model``). With millions of points the model ECDF at the observed
stars is counted exactly (``ks_model_ecdf = "count"``, ``stats.ecdf_at``); at 1M stars that agrees with the
interpolated ECDF to 7e-5 in KS (``logs/validate_count.out``).

For every (run, age) of the base grid below the run's age limit whose 32k minimum KS is at most FACTOR x
the grid's best: start at that minimum, evaluate the 8 neighbouring (mu, E) lattice cells, move to the best
while it lowers the KS, stop at a local minimum (as ks_refine/refine.py). Cells are replayed from
RandomState(crn_seed).

  refine4m.py OUT init BASE_REFINED_DIR AGE_LIMITS [--factor 2.5] [--per-unit 6]
        -> OUT/starts.csv (unit, run, age, mu0, e0, ks32k); OUT/setup.json must exist
  refine4m.py OUT worker [--max-hours H] [--max-steps 30]
  refine4m.py OUT status
  refine4m.py OUT merge  -> OUT/combined_2d_ks_results.csv (every cell evaluated, 4M sCMDs only)

Outputs: OUT/units/<unit>.csv (written after each (run, age)), OUT/claims/<unit>/ (owner, error).
"""
import argparse
import glob
import json
import os
import socket
import sys
import time
import traceback

import numpy as np
import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
sys.path.insert(0, f"{B}/m15_scmd4m/code")
COLS = ["simulation_id", "age_gyr", "distance_modulus", "extinction_ebv", "ks_statistic", "js_divergence",
        "mean_diff", "p95_diff", "n_observed", "n_simulated", "step"]


def init(a):
    lim = {int(k): v for k, v in json.load(open(a.age_limits)).items()}
    g = pd.read_csv(f"{a.base}/combined_2d_ks_results.csv",
                    usecols=["simulation_id", "age_gyr", "distance_modulus", "extinction_ebv", "ks_statistic"])
    g = g[g.age_gyr < g.simulation_id.map(lim) - 1e-6]
    b = g.loc[g.groupby(["simulation_id", "age_gyr"]).ks_statistic.idxmin()]
    b = b[b.ks_statistic <= a.factor * b.ks_statistic.min()].sort_values(["simulation_id", "age_gyr"])
    rows = []
    for run, r in b.groupby("simulation_id"):
        for k in range(0, len(r), a.per_unit):
            for _, x in r.iloc[k:k + a.per_unit].iterrows():
                rows.append({"unit": f"{run}_{k // a.per_unit}", "run": int(run), "age": float(x.age_gyr),
                             "mu0": float(x.distance_modulus), "e0": float(x.extinction_ebv),
                             "ks32k": float(x.ks_statistic)})
    s = pd.DataFrame(rows)
    s.to_csv(f"{a.out}/starts.csv", index=False)
    print(f"{len(s)} (run, age) starts from {s.run.nunique()} runs in {s.unit.nunique()} units "
          f"(32k KS <= {a.factor} x {b.ks_statistic.min():.5f})")


def lattice(S):
    dm, eb = np.round(np.arange(*S.dm_grid), 4), np.round(np.arange(*S.ebv_grid), 4)
    return dm, eb


def descend(ctx, make, dm, eb, i0, j0, max_steps):
    cache, rows = {}, []

    def ks(i, j, step):
        if (i, j) not in cache:
            cat = make(np.random.RandomState(ctx.s.crn_seed), dm[i], eb[j])
            if cat is None or len(cat) < 100:
                cache[(i, j)] = np.inf
            else:
                out = ctx.stage3(cat)
                cache[(i, j)] = out["ks_statistic"]
                rows.append({**{k: out[k] for k in ("ks_statistic", "js_divergence", "mean_diff", "p95_diff",
                                                    "n_observed", "n_simulated")},
                             "distance_modulus": dm[i], "extinction_ebv": eb[j], "step": step})
        return cache[(i, j)]

    cur, best, step = (i0, j0), ks(i0, j0, 0), 0
    while step < max_steps:
        step += 1
        nb = [(cur[0] + di, cur[1] + dj) for di in (-1, 0, 1) for dj in (-1, 0, 1)
              if (di, dj) != (0, 0) and 0 <= cur[0] + di < len(dm) and 0 <= cur[1] + dj < len(eb)]
        vals = [ks(i, j, step) for i, j in nb]
        k = int(np.argmin(vals))
        if vals[k] < best:
            cur, best = nb[k], vals[k]
        else:
            break
    return rows


def worker(a):
    from isochronetoolbox.scmd import Setup, context
    S = Setup(**{k: v for k, v in json.load(open(f"{a.out}/setup.json")).items() if k in Setup.__dataclass_fields__})
    dm, eb = lattice(S)
    starts = pd.read_csv(f"{a.out}/starts.csv")
    units = list(dict.fromkeys(starts.unit))
    rng = np.random.default_rng(int(os.environ.get("SLURM_ARRAY_TASK_ID", "0")))
    rng.shuffle(units)
    tag = f"{socket.gethostname()} {os.environ.get('SLURM_ARRAY_JOB_ID', '-')}_{os.environ.get('SLURM_ARRAY_TASK_ID', '-')}"
    t_end = time.time() + 3600 * a.max_hours
    os.makedirs(f"{a.out}/claims", exist_ok=True)
    os.makedirs(f"{a.out}/units", exist_ok=True)
    ctx = None
    for u in units:
        if time.time() > t_end:
            break
        if os.path.exists(f"{a.out}/units/{u}.csv"):
            continue
        try:
            os.mkdir(f"{a.out}/claims/{u}")
        except FileExistsError:
            continue
        open(f"{a.out}/claims/{u}/owner", "w").write(f"{tag} {time.strftime('%F %T')}\n")
        try:
            if ctx is None:
                ctx = context(S)
            t0, out = time.time(), []
            for _, r in starts[starts.unit == u].iterrows():
                _, make = ctx.model(int(r.run), float(r.age))
                i0, j0 = int(np.argmin(np.abs(dm - r.mu0))), int(np.argmin(np.abs(eb - r.e0)))
                for row in descend(ctx, make, dm, eb, i0, j0, a.max_steps):
                    out.append({"simulation_id": int(r.run), "age_gyr": float(r.age), **row})
                pd.DataFrame(out, columns=COLS).to_csv(f"{a.out}/units/{u}.partial", index=False)
            os.replace(f"{a.out}/units/{u}.partial", f"{a.out}/units/{u}.csv")
            print(f"{u}: {len(out)} cells, {time.time() - t0:.0f} s", flush=True)
        except Exception:
            open(f"{a.out}/claims/{u}/error", "w").write(traceback.format_exc())
            print(f"{u}: error\n{traceback.format_exc()}", flush=True)


def status(a):
    units = set(pd.read_csv(f"{a.out}/starts.csv").unit)
    done = {os.path.basename(p)[:-4] for p in glob.glob(f"{a.out}/units/*.csv")}
    claims = {os.path.basename(p) for p in glob.glob(f"{a.out}/claims/*")}
    err = len(glob.glob(f"{a.out}/claims/*/error"))
    print({"done": len(done & units), "claimed": len((claims - done) & units),
           "todo": len(units - claims - done), "errors": err})


def merge(a):
    parts = [pd.read_csv(p) for p in sorted(glob.glob(f"{a.out}/units/*.csv"))]
    m = pd.concat(parts, ignore_index=True)
    m.to_csv(f"{a.out}/combined_2d_ks_results.csv", index=False)
    s = pd.read_csv(f"{a.out}/starts.csv")
    b = m.loc[m.groupby(["simulation_id", "age_gyr"]).ks_statistic.idxmin()]
    j = b.merge(s, left_on=["simulation_id", "age_gyr"], right_on=["run", "age"])
    moved = ((j.distance_modulus - j.mu0).abs() > 1e-6) | ((j.extinction_ebv - j.e0).abs() > 1e-6)
    print(f"{len(m)} cells evaluated for {len(b)} (run, age) of {len(s)} starts; cells per (run, age) median "
          f"{m.groupby(['simulation_id', 'age_gyr']).size().median():.0f}; minimum moved from the 32k one in "
          f"{moved.mean() * 100:.0f}%; KS 4M - 32k at the minimum: median {np.median(j.ks_statistic - j.ks32k):+.5f}")
    r = b.loc[b.ks_statistic.idxmin()]
    print(f"best: run {int(r.simulation_id)} {r.age_gyr} Gyr mu {r.distance_modulus} E {r.extinction_ebv} "
          f"KS {r.ks_statistic:.5f} (32k grid best {s.ks32k.min():.5f})")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("out")
    sp = p.add_subparsers(dest="cmd", required=True)
    q = sp.add_parser("init"); q.add_argument("base"); q.add_argument("age_limits")
    q.add_argument("--factor", type=float, default=2.5); q.add_argument("--per-unit", type=int, default=6)
    q = sp.add_parser("worker"); q.add_argument("--max-hours", type=float, default=18.0)
    q.add_argument("--max-steps", type=int, default=30)
    sp.add_parser("status"); sp.add_parser("merge")
    a = p.parse_args()
    a.out = os.path.abspath(a.out)
    {"init": init, "worker": worker, "status": status, "merge": merge}[a.cmd](a)
