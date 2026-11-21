"""Robustness of the expected-ECDF split null to the AST split, and AST-list size scaling.

usage: python split_check.py MODE   (MODE = seed1 | seed2 | quarter)
  seedS  : halves(seed=S): pool A' (12 sCMDs), pool B' (12), null sets A'' (12), B'' (12)
  quarter: two disjoint random quarters Q1, Q2 of the artificial stars, same layout
For each pair (X, Y): null sets of Y vs pool X (split), null sets of X vs pool X (sampling only),
and the data vs each pool.  Statistic: P95 |ECDF_own - ECDF_pool| (exact).
"""
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, '/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/expected-ecdf')
from ecdf_tools import AGE, OUT, RUN, BigECDF, make_points, setup_context  # noqa: E402
from isochronetoolbox.scmd.ast import ArtificialStars  # noqa: E402
from isochronetoolbox.scmd.stats import rank2d  # noqa: E402

mode = sys.argv[1]
t0 = time.time()
S, ctx = setup_context()
full = ctx.ast
if mode.startswith('seed'):
    X, Y = full.halves(seed=int(mode[4:]))
    base = 600000 + 100000 * int(mode[4:])
else:
    n = full.n_ast
    u = np.random.RandomState(11).random_sample(n)
    X = ArtificialStars(*full._args, rows=u < 0.25)
    Y = ArtificialStars(*full._args, rows=(u >= 0.25) & (u < 0.5))
    base = 900000
_, make = ctx.model(RUN, AGE)
obs, obs_cdf = ctx.obs3, ctx.obs3_cdf
N = len(obs)


def gen(ast, seeds):
    ctx.ast = ast
    out = [make_points(ctx, S, make, s) for s in seeds]
    ctx.ast = full
    return out


pools = {'X': gen(X, [base + k for k in range(12)]), 'Y': gen(Y, [base + 1000 + k for k in range(12)])}
nulls = {'X': gen(X, [base + 2000 + k for k in range(12)]), 'Y': gen(Y, [base + 3000 + k for k in range(12)])}
print(f'{mode}: n_ast X {X.n_ast} Y {Y.n_ast}; generated in {time.time() - t0:.1f} s', flush=True)
E = {k: BigECDF(np.vstack(v)) for k, v in pools.items()}
rows = []
for k in 'XY':
    d = np.abs(obs_cdf - E[k](obs))
    rows.append({'mode': mode, 'kind': 'data', 'set': 'data', 'pool': k, 'i': -1, 'p95': np.percentile(d, 95),
                 'n_ast_pool': (X if k == 'X' else Y).n_ast})
for sname, sets in nulls.items():
    for i, x in enumerate(sets):
        own = rank2d(x) / N
        for k in 'XY':
            d = np.abs(own - E[k](x))
            rows.append({'mode': mode, 'kind': 'split' if k != sname else 'same', 'set': sname, 'pool': k, 'i': i,
                         'p95': np.percentile(d, 95), 'n_ast_pool': (X if k == 'X' else Y).n_ast})
df = pd.DataFrame(rows)
df.to_csv(f'{OUT}/rows/split_check_{mode}.csv', index=False)
for (kind, pool), g in df.groupby(['kind', 'pool']):
    print(f'{mode} {kind:5s} pool {pool}: p95 {g.p95.mean():.5f} +- {g.p95.std(ddof=1) if len(g) > 1 else 0:.5f} (n={len(g)})',
          flush=True)
print(f'{mode}: done in {time.time() - t0:.1f} s', flush=True)
