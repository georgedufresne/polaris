"""One worker of the orientation audit.

roles
  ref   : the split-AST bootstrap reference (grid seed, ctx.ast = half_a), saved to ref_split.npy
  grid  : data vs the grid's own best-cell sCMD, regenerated with the GRID's code snapshot (full ASTs)
  data  : data vs fresh realisations RandomState(i) of the best cell (latest code, full ASTs)
  null  : reference vs realisations RandomState(i) with ctx.ast = half_b (latest code)
"""
import argparse
import json
import os
import sys
import time

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("role", choices=["ref", "grid", "data", "null"])
ap.add_argument("--ids", default="")
ap.add_argument("--out", required=True)
ap.add_argument("--save", type=int, default=0, help="save per-point arrays for the first N realisations")
ap.add_argument("--validate", type=int, default=1, help="literal checks on the first realisation")
a = ap.parse_args()

HERE = "/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/orientation"
G = "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02"
TB = G + "/code" if a.role == "grid" else "/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox"
sys.path.insert(0, TB)
sys.path.insert(0, HERE)
import isochronetoolbox  # noqa: E402
from isochronetoolbox.scmd import Setup, context  # noqa: E402
from isochronetoolbox.scmd.pipeline import _window, regenerate  # noqa: E402
from isochronetoolbox.scmd.stats import ks_2d, rank2d  # noqa: E402
import orientlib as L  # noqa: E402

RUN, AGE, MU, E = 10099, 11.4, 13.38, 0.015
logf = open(a.out + ".log", "w")


def P(*x):
    print(*x, file=logf, flush=True)


t0 = time.time()
S = Setup(**{k: v for k, v in json.load(open(G + "/setup.json")).items() if k in Setup.__dataclass_fields__})
ctx = context(S)
P("package", isochronetoolbox.__file__, f"context {time.time() - t0:.1f} s, n_obs {len(ctx.obs3)}")
ids = [int(x) for x in a.ids.split(",") if x] if "," in a.ids or a.ids.isdigit() else \
      (list(range(*[int(x) for x in a.ids.split(":")])) if ":" in a.ids else [])


def win(cat):
    pts = ctx.points(cat)
    return pts[_window(pts, S.ks_window)]


def emit(fh, rec):
    fh.write(json.dumps(rec) + "\n")
    fh.flush()


if a.role == "ref":
    t = time.time()
    half_a, half_b = ctx.ast.halves()
    P(f"halves {time.time() - t:.1f} s: {half_a.n_ast} + {half_b.n_ast} of {ctx.ast.n_ast}")
    ctx.ast = half_a
    t = time.time()
    from isochronetoolbox.scmd.bootstrap import reference
    ref = reference(ctx, RUN, AGE, MU, E)
    P(f"reference {time.time() - t:.1f} s, {len(ref)} stars")
    np.save(os.path.join(HERE, "ref_split.npy"), ref)
    sys.exit(0)

if a.role == "null":
    half_a, half_b = ctx.ast.halves()
    ctx.ast = half_b
    ref = np.load(os.path.join(HERE, "ref_split.npy"))
else:
    ref = ctx.obs3
t = time.time()
F_ref = L.ref_ecdfs(ref, rank2d)
P(f"reference ECDFs {time.time() - t:.1f} s; bb equals pipeline obs3_cdf: "
  f"{np.array_equal(F_ref['bb'], ctx.obs3_cdf) if a.role != 'null' else 'n/a'}")

make = None if a.role == "grid" else ctx.model(RUN, AGE)[1]   # the grid snapshot has no Context.model
with open(a.out, "w") as fh:
    if a.role == "grid":
        t = time.time()
        sim = win(regenerate(ctx, RUN, AGE, MU, E))
        tg = time.time() - t
        ks_pipe = float(ks_2d(ctx.obs3, sim, ctx.obs3_cdf)[0])
        rec, arr = L.score(ref, F_ref, sim, rank2d, validate=bool(a.validate))
        rec.update(role="grid", id=-1, ks_pipeline=ks_pipe, t_gen=tg, t_total=time.time() - t)
        emit(fh, rec)
        if a.save:
            np.savez_compressed(a.out + ".grid.npz", ref=ref, sim=sim, **arr)
        P(f"grid sCMD: pipeline ks_2d {ks_pipe:.10f}, bb_p95 {rec['bb_p95']:.10f}, {time.time() - t:.1f} s")
    else:
        for j, i in enumerate(ids):
            t = time.time()
            sim = win(make(np.random.RandomState(i), MU, E))
            tg = time.time() - t
            val = bool(a.validate) and j == 0
            rec, arr = L.score(ref, F_ref, sim, rank2d, validate=val)
            if val:
                rec["ks_pipeline"] = float(ks_2d(ref, sim, F_ref["bb"])[0])
            rec.update(role=a.role, id=i, t_gen=tg, t_total=time.time() - t)
            emit(fh, rec)
            if j < a.save:
                np.savez_compressed(f"{a.out}.{i}.npz", ref=ref, sim=sim, **arr)
            P(f"{a.role} {i}: bb_p95 {rec['bb_p95']:.6f} bf {rec['bf_p95']:.6f} rb {rec['rb_p95']:.6f} "
              f"rf {rec['rf_p95']:.6f}  {time.time() - t:.1f} s")
P(f"done {time.time() - t0:.1f} s")
