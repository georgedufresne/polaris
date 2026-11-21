"""Does the split null's level depend on which random AST split is used?  halves(seed=S):
4 references (half_a', seeds 910001-4) x 12 realisations (half_b', seeds 0-11)."""
import sys
import numpy as np
sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/reference-selection")
from common import AGE, EVI, MU, OUT, RUN, Log, ks_many, load_ctx, prep_target, setup_path, window_points
setup_path("latest")
SEED = int(sys.argv[1])
log = Log(f"split_seed{SEED}")
S, ctx = load_ctx()
full = ctx.ast
ha, hb = full.halves(seed=SEED)
_, make = ctx.model(RUN, AGE)
ctx.ast = ha
refs = [prep_target(window_points(ctx, S, make(np.random.RandomState(910000 + k), MU, EVI))) for k in range(1, 5)]
ctx.ast = hb
Y = []
for j in range(12):
    r = ks_many(window_points(ctx, S, make(np.random.RandomState(j), MU, EVI)), refs)
    Y.append([v["ks"] for v in r])
Y = np.array(Y)
np.save(f"{OUT}/split_seed{SEED}_ks.npy", Y)
log(f"split seed {SEED}: per-reference mean " + " ".join(f"{m:.5f}" for m in Y.mean(0)) +
    f"; per-reference sd " + " ".join(f"{s:.5f}" for s in Y.std(0, ddof=1)) + f"; grand mean {Y.mean():.5f}")
