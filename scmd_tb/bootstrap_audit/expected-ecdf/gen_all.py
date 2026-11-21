"""Generate the point sets (colour, F606W) in the KS window at the best cell.

usage: python gen_all.py GROUP
  bigA  : 24 realisations, AST half a (seeds 100000+k)  -> pooled 'expected' model A
  bigB  : 24 realisations, AST half b (seeds 400000+k)  -> pooled model B (symmetry check)
  bigF  : 24 realisations, full AST list (seeds 300000+k) -> pooled model F (as the grid)
  null  : 48 realisations, AST half b (seeds 0..47)  [null 'data' stand-ins, as the bootstrap]
          24 realisations, AST half a (seeds 200000+k) [sampling-only control]
"""
import os
import sys
import time

import numpy as np

sys.path.insert(0, '/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/expected-ecdf')
from ecdf_tools import AGE, OUT, RUN, make_points, setup_context  # noqa: E402

group = sys.argv[1]
os.makedirs(f'{OUT}/sets', exist_ok=True)
t0 = time.time()
S, ctx = setup_context()
full = ctx.ast
half_a, half_b = full.halves()
_, make = ctx.model(RUN, AGE)
plan = {
    'bigA': [('bigA', half_a, [100000 + k for k in range(24)])],
    'bigB': [('bigB', half_b, [400000 + k for k in range(24)])],
    'bigF': [('bigF', full, [300000 + k for k in range(24)])],
    'null': [('nullB', half_b, list(range(48))), ('nullA', half_a, [200000 + k for k in range(24)])],
}[group]
for name, ast, seeds in plan:
    ctx.ast = ast
    pts = [make_points(ctx, S, make, s) for s in seeds]
    ctx.ast = full
    n = np.array([len(p) for p in pts])
    np.savez(f'{OUT}/sets/{name}.npz', pts=np.stack(pts) if len(set(n)) == 1 else np.array(pts, dtype=object),
             seeds=np.array(seeds), n=n, n_ast=ast.n_ast)
    print(f'{name}: {len(seeds)} sets, sizes {n.min()}-{n.max()}, n_ast {ast.n_ast}, '
          f'{time.time() - t0:.1f} s', flush=True)
