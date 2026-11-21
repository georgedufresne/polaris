"""Score data and null realisations against the pooled 'expected' models and N-vs-N.

usage: python score.py TASK NTASKS
Writes rows to OUT/rows/rows_TASK.csv.  Statistic definitions (all on |ECDF difference|):
  expected_exact : |ECDF_x(x) - ECDF_big(x)| at x's own points, ECDF_big exact (BigECDF)
  NvN_exact      : |ECDF_ref(ref) - ECDF_sim(ref)|, ECDF_sim exact cross count (/N_sim)
  NvN_interp     : the pipeline's ks_2d(ref, sim, ref_cdf) (LinearNDInterpolator, 0 outside hull)
'ref' is the set in the observed role (the data, or the bootstrap reference).
"""
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, '/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/expected-ecdf')
from ecdf_tools import OUT, BigECDF, setup_context  # noqa: E402
from isochronetoolbox.scmd.stats import ks_2d, rank2d  # noqa: E402

task, ntasks = int(sys.argv[1]), int(sys.argv[2])
os.makedirs(f'{OUT}/rows', exist_ok=True)
t0 = time.time()


def log(msg):
    print(f'[{time.time() - t0:7.1f} s] {msg}', flush=True)


def stats(own, other):
    d = np.abs(own - other)
    m = (own != 0) & (other != 0)
    dm = d[m]
    k = int(len(dm) * 0.99)
    return {'p95': float(np.percentile(d, 95)), 'p99': float(np.percentile(d, 99)),
            'max': float(d.max()), 'mean': float(d.mean()),
            'p99_ying': float(np.partition(dm, k)[k]), 'mean_signed': float(np.mean(own - other))}


S, ctx = setup_context()
obs, obs_cdf = ctx.obs3, ctx.obs3_cdf
N = len(obs)
sets = {k: np.load(f'{OUT}/sets/{k}.npz', allow_pickle=True) for k in ['bigA', 'bigB', 'bigF', 'nullB', 'nullA']}
P = {k: v['pts'] for k, v in sets.items()}
seeds = {k: v['seeds'] for k, v in sets.items()}
big = {'bigA': BigECDF(P['bigA'].reshape(-1, 2)), 'bigB': BigECDF(P['bigB'].reshape(-1, 2)),
       'bigF': BigECDF(P['bigF'].reshape(-1, 2)),
       'bigA_h1': BigECDF(P['bigA'][:12].reshape(-1, 2)), 'bigA_h2': BigECDF(P['bigA'][12:].reshape(-1, 2))}
log(f'big models built: ' + ', '.join(f'{k} {v.n}' for k, v in big.items()))
rows = []

# ---------------------------------------------------------------- data, expected version
if task == 0:
    variants = {'rank2d': (obs, obs_cdf)}
    Eo = BigECDF(obs)
    strict = Eo.count(obs)
    nonstrict = Eo.count(obs + 0.0004) - 1
    variants['midpoint'] = (obs, (strict + nonstrict) / 2 / N)
    log(f'obs grid check: max |1000 x - round| = {np.max(np.abs(obs * 1000 - np.round(obs * 1000))):.2e}; '
        f'rank2d - strict mean {np.mean(obs_cdf * N - strict):.1f}, nonstrict - strict mean {np.mean(nonstrict - strict):.1f}')
    for j in range(3):
        rs = np.random.RandomState(777 + j)
        dv, di = rs.uniform(-5e-4, 5e-4, N), rs.uniform(-5e-4, 5e-4, N)
        oj = np.column_stack([obs[:, 0] + dv - di, obs[:, 1] + dv])
        variants[f'jitter{j}'] = (oj, rank2d(oj) / N)
    for vname, (q, own) in variants.items():
        for bname, E in big.items():
            r = {'role': 'data', 'idx': -1, 'seed': -1, 'method': 'expected_exact', 'versus': bname,
                 'ref': 'data', 'obs_ecdf': vname}
            r.update(stats(own, E(q)))
            rows.append(r)
    # pools of increasing size (convergence of the expected model)
    for m in (1, 3, 6, 12, 24):
        E = BigECDF(P['bigA'][:m].reshape(-1, 2))
        r = {'role': 'data', 'idx': m, 'seed': -1, 'method': 'expected_exact', 'versus': f'bigA_first{m}',
             'ref': 'data', 'obs_ecdf': 'rank2d'}
        r.update(stats(obs_cdf, E(obs)))
        rows.append(r)
    # save the per-point difference against bigA for later diagnostics
    np.savez(f'{OUT}/data_vs_bigA.npz', obs=obs, obs_cdf=obs_cdf, big_cdf=big['bigA'](obs))
    log('data vs pooled models done')
    # data vs every individual member, exact cross count (N vs N)
    for bname in ('bigA', 'bigB', 'bigF'):
        for k in range(P[bname].shape[0]):
            Em = BigECDF(P[bname][k])
            r = {'role': 'data', 'idx': k, 'seed': int(seeds[bname][k]), 'method': 'NvN_exact',
                 'versus': f'{bname}_member', 'ref': 'data', 'obs_ecdf': 'rank2d'}
            r.update(stats(obs_cdf, Em(obs)))
            rows.append(r)
    log('data vs members (exact) done')

# data vs members, interpolated (the pipeline's ks_2d), split over tasks
for bname in ('bigA', 'bigF'):
    for k in range(task, P[bname].shape[0], ntasks):
        ks, d = ks_2d(obs, P[bname][k], obs_cdf)
        k99 = int(len(d) * 0.99)
        rows.append({'role': 'data', 'idx': k, 'seed': int(seeds[bname][k]), 'method': 'NvN_interp',
                     'versus': f'{bname}_member', 'ref': 'data', 'obs_ecdf': 'rank2d', 'p95': float(ks),
                     'p99': float(np.percentile(d, 99)), 'max': float(d.max()), 'mean': float(d.mean())})
log('data vs members (interp) done')

# ---------------------------------------------------------------- null realisations
refs_A = [(f'bigA_m{k}', P['bigA'][k]) for k in range(4)]    # half a references (N stars)
refs_B = [(f'bigB_m{k}', P['bigB'][k]) for k in range(4)]    # half b references
ref_cdf = {name: rank2d(p) / N for name, p in refs_A + refs_B}
for role, others, interp_refs in (('nullB', refs_A, ['bigA_m0', 'bigA_m1']),
                                  ('nullA', refs_B + refs_A, ['bigB_m0', 'bigA_m0'])):
    for i in range(task, P[role].shape[0], ntasks):
        x = P[role][i]
        own = rank2d(x) / N
        for bname, E in big.items():
            r = {'role': role, 'idx': i, 'seed': int(seeds[role][i]), 'method': 'expected_exact',
                 'versus': bname, 'ref': role, 'obs_ecdf': 'rank2d'}
            r.update(stats(own, E(x)))
            rows.append(r)
        Ex = BigECDF(x)
        for name, p in others:
            r = {'role': role, 'idx': i, 'seed': int(seeds[role][i]), 'method': 'NvN_exact',
                 'versus': role, 'ref': name, 'obs_ecdf': 'rank2d'}
            r.update(stats(ref_cdf[name], Ex(p)))
            rows.append(r)
        for name in interp_refs:
            p = dict(others)[name]
            ks, d = ks_2d(p, x, ref_cdf[name])
            rows.append({'role': role, 'idx': i, 'seed': int(seeds[role][i]), 'method': 'NvN_interp',
                         'versus': role, 'ref': name, 'obs_ecdf': 'rank2d', 'p95': float(ks),
                         'p99': float(np.percentile(d, 99)), 'max': float(d.max()), 'mean': float(d.mean())})
        log(f'{role} {i} done')

pd.DataFrame(rows).to_csv(f'{OUT}/rows/rows_{task:02d}.csv', index=False)
log(f'wrote {len(rows)} rows')
