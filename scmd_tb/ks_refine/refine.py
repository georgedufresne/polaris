#!/usr/bin/env python3
"""KS refinement of a finished sCMD grid: local descent on the 2D KS over the (mu, E) grid of every
(run, age), starting from the cells stage 2 passed to stage 3.

Why: stage 2 (fixed-bandwidth KDE JS) and stage 3 (2D KS) disagree on where the best (mu, E) is. In
m15_final, run 10392 at 13.6 Gyr, stage 2 kept 12 cells with KS 0.050-0.084 while (mu 15.48, E 0.060)
has KS 0.0150; that cell ranked 47th of 465 on the KDE JS.

For each (run, age) below the run's age limit: start from the previous age's optimum in the same unit
(and from the lowest-KS kept cell when that is already lower, or at the unit's first age), evaluate the 8 neighbouring grid cells, move to the best one while it
lowers the KS, stop at a local minimum. Catalogues are replayed exactly (crn: every cell is drawn from
RandomState(crn_seed)), and scored with the grid's own stage 3, so a refined row is identical to the row
the grid would have written for that cell.

  refine.py OUT worker [--age-limits JSON] [--max-hours H] [--first RUNS_FILE]
  refine.py OUT status
  refine.py OUT merge        -> OUT/combined_2d_ks_results_refined.csv (grid rows + refined rows)

Outputs: OUT/refine/units/<run>/<unit>.csv (every cell evaluated), OUT/refine/claims/, OUT/refine/logs/.
Optional OUT/refine/setup_override.json replaces setup fields for the descent (m15_final: ebv_grid
extended from 0.050 down to 0.000, the E(B-V) prior lowered on 2026-11-05).
"""
import argparse
import glob
import json
import os
import socket
import sys
import time

import numpy as np
import pandas as pd

OUT = os.path.abspath(sys.argv[1])
sys.path.insert(0, os.path.join(OUT, "code"))
from isochronetoolbox.scmd import Setup, context                                     # noqa: E402
from isochronetoolbox.scmd.pipeline import KS_COLUMNS, label, unit_key, units        # noqa: E402

R = os.path.join(OUT, "refine")
NBR = [(di, dj) for di in (-1, 0, 1) for dj in (-1, 0, 1) if (di, dj) != (0, 0)]


def load_setup():
    """The grid's setup, with OUT/refine/setup_override.json applied if present (e.g. a wider ebv_grid:
    the descent may then walk into cells the grid never had; the stage-2 cells still seed it)."""
    base = json.load(open(os.path.join(OUT, "setup.json")))
    ov = os.path.join(R, "setup_override.json")
    if os.path.exists(ov):
        base.update(json.load(open(ov)))
    return Setup(**{k: v for k, v in base.items() if k in Setup.__dataclass_fields__})


def all_units(s, first=None):
    runs = [int(x) for x in open(os.path.join(OUT, "runs.txt")).read().split()]
    if first:
        f = [int(x) for x in open(first).read().split()]
        runs = [r for r in f if r in runs] + [r for r in runs if r not in f]
    return units(runs, list(s.ages), s.chunk)


def out_path(run, key):
    return os.path.join(R, "units", str(run), key + ".csv")


def refine_unit(ctx, run, ages, limit, log, only=None):
    s = ctx.s
    key = unit_key(run, ages)
    ks = pd.read_csv(os.path.join(OUT, "units", str(run), key + ".ks.csv"))
    js = pd.read_csv(os.path.join(OUT, "units", str(run), key + ".js.csv"))
    rows, prev = [], None
    for age in ages:
        if (limit is not None and age >= limit - 1e-6) or (only is not None and age not in only):
            continue
        k_a = ks[(ks.age_gyr - age).abs() < 1e-6]
        j_a = js[(js.age_gyr - age).abs() < 1e-6].copy()
        if k_a.empty:
            log(f"{run} {age}: no stage-3 rows, skipped")
            continue
        j_a["rank"] = j_a.js_distance.rank(method="first").astype(int)
        jmap = {(round(d, 3), round(e, 3)): (v, r) for d, e, v, r in
                zip(j_a.distance_modulus, j_a.extinction_ebv, j_a.js_distance, j_a["rank"])}
        idx = lambda d, e: (int(np.argmin(np.abs(s.dm - d))), int(np.argmin(np.abs(s.ebv - e))))
        score = {idx(d, e): v for d, e, v in zip(k_a.distance_modulus, k_a.extinction_ebv, k_a.ks_statistic)}
        age_yr, make = ctx.model(run, age)
        n_eval = 0

        def ev(c):
            nonlocal n_eval
            if c in score:
                return score[c]
            i, j = c
            cat = make(np.random.RandomState(s.crn_seed), s.dm[i], s.ebv[j])
            n_eval += 1
            if cat is None:
                score[c] = np.inf
                return np.inf
            row = ctx.stage3(cat)
            name, a, d, e = label(age_yr, s.dm[i], s.ebv[j], s.name)
            jv, jr = jmap.get((round(d, 3), round(e, 3)), (np.nan, -1))
            row.update(simulation_id=run, filename=name, age_gyr=a, distance_modulus=d, extinction_ebv=e,
                       file_path="", plot_path="", kde_js=jv, kde_rank=jr)
            rows.append(row)
            score[c] = row["ks_statistic"]
            return score[c]

        start = min(score, key=score.get)
        ks0 = score[start]
        # the previous age's optimum is usually a few steps from this age's; the stage-2 cells are
        # descended from as well when they are already better than it (and always at a unit's first age)
        if prev is None:
            starts = [start]
        else:
            starts = [prev] + ([start] if ks0 < ev(prev) else [])
        for c in starts:
            cur = c
            ev(cur)
            while True:
                nb = [(cur[0] + di, cur[1] + dj) for di, dj in NBR
                      if 0 <= cur[0] + di < len(s.dm) and 0 <= cur[1] + dj < len(s.ebv)]
                best = min(nb, key=ev)
                if score[best] < score[cur]:
                    cur = best
                else:
                    break
        prev = min(score, key=score.get)
        log(f"{run} {age}: kept best {ks0:.5f} -> {score[prev]:.5f} at mu {s.dm[prev[0]]:.3f} "
            f"E {s.ebv[prev[1]]:.3f}; {n_eval} cells scored")
    return pd.DataFrame(rows, columns=KS_COLUMNS)


def worker(args):
    s = load_setup()
    assert s.crn, "refinement replays cells from crn_seed; this grid was not run with crn"
    lim = None
    if args.age_limits:
        lim = {int(k): v for k, v in json.load(open(args.age_limits)).items()}
    ctx = context(s)
    t_end = time.time() + 3600 * args.max_hours
    tag = f"{socket.gethostname()}:{os.environ.get('SLURM_ARRAY_JOB_ID', '')}_{os.environ.get('SLURM_ARRAY_TASK_ID', '')}"
    log = lambda m: print(m, flush=True)
    for run, ages in all_units(s, args.first):
        if time.time() > t_end:
            log("max hours reached, stopping")
            break
        key = unit_key(run, ages)
        if os.path.exists(out_path(run, key)):
            continue
        try:
            os.mkdir(os.path.join(R, "claims", key))
        except FileExistsError:
            continue
        open(os.path.join(R, "claims", key, "owner"), "w").write(f"{tag} {time.strftime('%F %T')}\n")
        t0 = time.time()
        try:
            df = refine_unit(ctx, run, ages, lim.get(run) if lim else None, log)
        except Exception:
            import traceback
            open(os.path.join(R, "claims", key, "error"), "w").write(traceback.format_exc())
            log(f"{key}: FAILED\n{traceback.format_exc()}")
            continue
        os.makedirs(os.path.dirname(out_path(run, key)), exist_ok=True)
        tmp = out_path(run, key) + f".tmp{os.getpid()}"
        df.to_csv(tmp, index=False)
        os.replace(tmp, out_path(run, key))
        log(f"{key}: {len(df)} rows, {time.time() - t0:.0f} s")


def status(args):
    s = load_setup()
    us = all_units(s)
    done = sum(os.path.exists(out_path(r, unit_key(r, a))) for r, a in us)
    claims = {os.path.basename(p) for p in glob.glob(os.path.join(R, "claims", "*"))}
    err = len(glob.glob(os.path.join(R, "claims", "*", "error")))
    claimed = sum(unit_key(r, a) in claims and not os.path.exists(out_path(r, unit_key(r, a))) for r, a in us)
    print({"done": done, "claimed": claimed, "todo": len(us) - done - claimed, "errors": err})


def merge(args):
    s = load_setup()
    us = all_units(s)
    miss = [unit_key(r, a) for r, a in us if not os.path.exists(out_path(r, unit_key(r, a)))]
    if miss:
        sys.exit(f"{len(miss)} units not refined yet, e.g. {miss[:3]}")
    g = pd.read_csv(os.path.join(OUT, "combined_2d_ks_results.csv")).assign(source="grid")
    parts = [pd.read_csv(out_path(r, unit_key(r, a))) for r, a in us]
    f = pd.concat([p for p in parts if len(p)], ignore_index=True).assign(source="refine")
    m = pd.concat([g, f], ignore_index=True)
    k = ["simulation_id", "age_gyr", "distance_modulus", "extinction_ebv"]
    m = m.drop_duplicates(subset=k, keep="first")
    m.to_csv(os.path.join(OUT, "combined_2d_ks_results_refined.csv"), index=False)
    gb = g.groupby(["simulation_id", "age_gyr"]).ks_statistic.min()
    mb = m.groupby(["simulation_id", "age_gyr"]).ks_statistic.min().reindex(gb.index)
    print(f"grid rows {len(g)}, refined rows added {len(m) - len(g)}; (run, age) cells improved: "
          f"{(mb < gb - 1e-12).sum()} of {len(gb)}; median improvement {np.median(gb - mb):.5f}")
    print(f"best KS: grid {g.ks_statistic.min():.5f} -> refined {m.ks_statistic.min():.5f}")
    print(os.path.join(OUT, "combined_2d_ks_results_refined.csv"))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("out")
    p.add_argument("command", choices=["worker", "status", "merge"])
    p.add_argument("--age-limits")
    p.add_argument("--max-hours", type=float, default=10.0)
    p.add_argument("--first", help="file of runs to refine first (validation runs)")
    a = p.parse_args()
    for d in ("units", "claims", "logs"):
        os.makedirs(os.path.join(R, d), exist_ok=True)
    {"worker": worker, "status": status, "merge": merge}[a.command](a)
