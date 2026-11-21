"""Phase 1: build the reference sCMDs (KS-window points) and save them as .npy.

usage: phase1_refs.py official_split | official_shared | fresh
  official_split : R0 = bootstrap.reference with ctx.ast = half_a (grid seed), as the
                   published split null
  official_shared: S0 = bootstrap.reference with the full AST list (grid seed), as the
                   published shared null (latest code: ridge through the first n_trace stars)
  fresh          : R1..R9 (half_a, seeds 900001..900009), S1..S5 (full, seeds 900101..900105)
"""
import sys

import numpy as np

sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/reference-selection")
from common import AGE, EVI, MU, OUT, RUN, Log, load_ctx, setup_path, window_points  # noqa: E402

setup_path("latest")
from isochronetoolbox.scmd.bootstrap import reference  # noqa: E402

mode = sys.argv[1]
log = Log("phase1_" + mode)
S, ctx = load_ctx()
log("context ready; n_obs3", len(ctx.obs3), "n_trace", ctx.rcut.n_trace, "n_ast", ctx.ast.n_ast)
full = ctx.ast
if mode in ("official_split", "fresh"):
    half_a, half_b = full.halves()
    log("halves", half_a.n_ast, half_b.n_ast)

if mode == "official_split":
    ctx.ast = half_a
    ref = reference(ctx, RUN, AGE, MU, EVI)
    np.save(f"{OUT}/ref_R0.npy", ref)
    log("R0 (official split reference, half_a, grid seed):", len(ref), "points")
elif mode == "official_shared":
    ctx.ast = full
    ref = reference(ctx, RUN, AGE, MU, EVI)
    np.save(f"{OUT}/ref_S0.npy", ref)
    log("S0 (official shared reference, full AST, grid seed):", len(ref), "points")
elif mode == "fresh":
    _, make = ctx.model(RUN, AGE)
    ctx.ast = half_a
    for k in range(1, 10):
        p = window_points(ctx, S, make(np.random.RandomState(900000 + k), MU, EVI))
        np.save(f"{OUT}/ref_R{k}.npy", p)
        log(f"R{k} (half_a, seed {900000 + k}):", len(p))
    ctx.ast = full
    for k in range(1, 6):
        p = window_points(ctx, S, make(np.random.RandomState(900100 + k), MU, EVI))
        np.save(f"{OUT}/ref_S{k}.npy", p)
        log(f"S{k} (full, seed {900100 + k}):", len(p))
log("done")
