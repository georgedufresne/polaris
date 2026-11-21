"""Summaries of the expected-ECDF audit rows."""
import glob

import numpy as np
import pandas as pd

OUT = '/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/expected-ecdf'
df = pd.concat([pd.read_csv(p) for p in sorted(glob.glob(f'{OUT}/rows/rows_*.csv'))], ignore_index=True)
df.to_csv(f'{OUT}/all_rows.csv', index=False)
pd.set_option('display.width', 200, 'display.max_rows', 200)


def ms(x):
    x = np.asarray(x, float)
    return f'{x.mean():.5f} +- {x.std(ddof=1):.5f} (n={len(x)}, ddof=1; 3sig {x.mean() + 3 * x.std(ddof=0):.5f} with ddof=0)'


print('=== DATA vs pooled models (expected_exact) ===')
d = df[(df.role == 'data') & (df.method == 'expected_exact')]
print(d[['versus', 'obs_ecdf', 'p95', 'p99', 'max', 'mean', 'p99_ying', 'mean_signed']].to_string(index=False))

print('\n=== DATA vs individual N-star models (N vs N) ===')
for (m, v), g in df[(df.role == 'data') & (df.method.str.startswith('NvN'))].groupby(['method', 'versus']):
    print(f'{m:11s} {v:12s} p95 {ms(g.p95)}   max {g["max"].mean():.4f}   p99 {g.p99.mean():.4f}')

print('\n=== NULL realisations vs pooled models (expected_exact) ===')
for (role, v), g in df[(df.role != 'data') & (df.method == 'expected_exact')].groupby(['role', 'versus']):
    print(f'{role:6s} vs {v:8s} p95 {ms(g.p95)}  | p99 {g.p99.mean():.5f} +- {g.p99.std(ddof=1):.5f} '
          f'| max {g["max"].mean():.4f} +- {g["max"].std(ddof=1):.4f} | p99_ying {g.p99_ying.mean():.5f} +- {g.p99_ying.std(ddof=1):.5f}')

print('\n=== NULL N vs N (reference in the observed role) ===')
for (role, m, ref), g in df[(df.role != 'data') & (df.method.str.startswith('NvN'))].groupby(['role', 'method', 'ref']):
    print(f'{role:6s} {m:11s} ref {ref:8s} p95 {ms(g.p95)}  | max {g["max"].mean():.4f} +- {g["max"].std(ddof=1):.4f}')
