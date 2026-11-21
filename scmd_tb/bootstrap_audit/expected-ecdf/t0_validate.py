"""Timing and validation: context, one sCMD, rank2d, and BigECDF against rank2d with a mask."""
import sys
import time

import numpy as np

sys.path.insert(0, '/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/expected-ecdf')
from ecdf_tools import AGE, RUN, BigECDF, make_points, setup_context  # noqa: E402
from isochronetoolbox.scmd.stats import rank2d  # noqa: E402

t = time.time()
S, ctx = setup_context()
print(f'context {time.time() - t:.1f} s; obs3 {ctx.obs3.shape}; n_target {ctx.n_target}; '
      f'n_ast {ctx.ast.n_ast}; n_trace {ctx.rcut.n_trace}', flush=True)
full = ctx.ast
t = time.time()
half_a, half_b = full.halves()
print(f'halves {time.time() - t:.1f} s; n_ast a {half_a.n_ast} b {half_b.n_ast}', flush=True)
_, make = ctx.model(RUN, AGE)
ctx.ast = half_a
t = time.time()
P = [make_points(ctx, S, make, 900000 + k) for k in range(3)]
print(f'3 sCMDs (half a) {time.time() - t:.1f} s; sizes {[len(p) for p in P]}', flush=True)
ctx.ast = full

t = time.time()
r = rank2d(P[0])
print(f'rank2d 47k: {time.time() - t:.2f} s', flush=True)

big = np.vstack(P)
q = ctx.obs3
t = time.time()
E = BigECDF(big)
t1 = time.time()
c = E.count(q)
print(f'BigECDF build {t1 - t:.2f} s, query {time.time() - t1:.2f} s (big {len(big)}, q {len(q)})', flush=True)

t = time.time()
allp = np.vstack([big, q])
mask = np.zeros(len(allp), bool)
mask[:len(big)] = True
c_ref = rank2d(allp, mask)[~mask]
print(f'rank2d masked on {len(allp)}: {time.time() - t:.1f} s', flush=True)
print('max |BigECDF - rank2d masked| (counts):', int(np.max(np.abs(c - c_ref))),
      ' n differing:', int(np.sum(c != c_ref)), flush=True)

# brute force on a random subset of queries
rs = np.random.RandomState(1)
sub = rs.choice(len(q), 300, replace=False)
bf = np.array([np.sum((big[:, 0] < q[j, 0]) & (big[:, 1] < q[j, 1])) for j in sub])
print('max |BigECDF - brute force| on 300 queries:', int(np.max(np.abs(c[sub] - bf))), flush=True)
# the realisation's own points as queries (ties with themselves): exclude-self check
own = P[0]
c_own = E.count(own)
bf_own = np.array([np.sum((big[:, 0] < own[j, 0]) & (big[:, 1] < own[j, 1])) for j in sub])
print('own-points queries, max |BigECDF - brute|:', int(np.max(np.abs(c_own[sub] - bf_own))), flush=True)
# rank2d of a single set vs brute force strict count
bfr = np.array([np.sum((own[:, 0] < own[j, 0]) & (own[:, 1] < own[j, 1])) for j in sub])
print('rank2d(own) vs strict brute force, max diff:', int(np.max(np.abs(r[sub] - bfr))), flush=True)
# observed ECDF from rank2d vs strict brute force (ties in the rounded data)
o = ctx.obs3
bfo = np.array([np.sum((o[:, 0] < o[j, 0]) & (o[:, 1] < o[j, 1])) for j in sub])
bfo_le = np.array([np.sum((o[:, 0] <= o[j, 0]) & (o[:, 1] <= o[j, 1])) - 1 for j in sub])
ro = np.rint(ctx.obs3_cdf * len(o)).astype(int)
print('obs rank2d vs strict brute: mean diff', float(np.mean(ro[sub] - bfo)), 'max', int(np.max(np.abs(ro[sub] - bfo))),
      '| vs non-strict: mean diff', float(np.mean(ro[sub] - bfo_le)), flush=True)
print('distinct obs colours', len(np.unique(o[:, 0])), 'distinct obs mags', len(np.unique(o[:, 1])))
