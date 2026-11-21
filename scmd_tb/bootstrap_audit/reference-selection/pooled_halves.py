"""Why does the split null depend so much on the reference?  Pool 9 sCMDs per AST list
(~427k stars, sampling noise ~3x lower) and compare the pools with the same KS:
a systematic difference between half_a and half_b sCMDs survives pooling, noise does not.
Also data vs a pooled full-list model (model noise ~3x lower)."""
import sys
import time

import numpy as np

sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/reference-selection")
from common import AGE, EVI, MU, OUT, RUN, Log, ks_many, load_ctx, prep_target, setup_path, window_points  # noqa: E402

setup_path("latest")
log = Log("pooled_halves")
S, ctx = load_ctx()
full = ctx.ast
half_a, half_b = full.halves()
_, make = ctx.model(RUN, AGE)


def pool(ast, seeds):
    ctx.ast = ast
    return np.vstack([window_points(ctx, S, make(np.random.RandomState(s), MU, EVI)) for s in seeds])


P = {
    "A1 (half_a, R1-R9)": np.vstack([np.load(f"{OUT}/ref_R{k}.npy") for k in range(1, 10)]),
    "A2 (half_a, seeds 700001-9)": pool(half_a, range(700001, 700010)),
    "B1 (half_b, seeds 0-8)": pool(half_b, range(0, 9)),
    "B2 (half_b, seeds 9-17)": pool(half_b, range(9, 18)),
    "F1 (full, S1-S5 + 900106-9)": np.vstack([np.load(f"{OUT}/ref_S{k}.npy") for k in range(1, 6)]
                                             + [pool(full, range(900106, 900110))]),
    "F2 (full, seeds 900201-9)": pool(full, range(900201, 900210)),
}
ctx.ast = full
log("pools built", {k: len(v) for k, v in P.items()})
T = {k: prep_target(v) for k, v in P.items()}
log("target ECDFs done")
T["data (47k)"] = prep_target(ctx.obs3, ctx.obs3_cdf)
pairs = [("A1 (half_a, R1-R9)", "A2 (half_a, seeds 700001-9)"),
         ("B1 (half_b, seeds 0-8)", "B2 (half_b, seeds 9-17)"),
         ("F1 (full, S1-S5 + 900106-9)", "F2 (full, seeds 900201-9)"),
         ("A1 (half_a, R1-R9)", "B1 (half_b, seeds 0-8)"),
         ("A2 (half_a, seeds 700001-9)", "B2 (half_b, seeds 9-17)"),
         ("A1 (half_a, R1-R9)", "F1 (full, S1-S5 + 900106-9)"),
         ("B1 (half_b, seeds 0-8)", "F2 (full, seeds 900201-9)"),
         ("data (47k)", "F1 (full, S1-S5 + 900106-9)"),
         ("data (47k)", "F2 (full, seeds 900201-9)")]
# group by the interpolated (sim) set so each Delaunay is built once
by_sim = {}
for ref, sim in pairs:
    by_sim.setdefault(sim, []).append(ref)
for sim, refs in by_sim.items():
    res = ks_many(P[sim], [T[r] for r in refs])
    for r, v in zip(refs, res):
        log(f"KS(ref = {r}, sim = {sim}) = {v['ks']:.6f}   p99 {v['p99']:.6f}  mean|d| {v['mean']:.6f}")
log("done")
