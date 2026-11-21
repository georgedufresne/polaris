"""Separations, pass/fail and the variance decomposition from all_rows.csv."""
import json

import numpy as np
import pandas as pd

OUT = '/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/expected-ecdf'
df = pd.read_csv(f'{OUT}/all_rows.csv')


def null(role, method, versus=None, ref=None, col='p95'):
    g = df[(df.role == role) & (df.method == method)]
    if versus is not None:
        g = g[g.versus == versus]
    if ref is not None:
        g = g[g.ref == ref]
    x = g[col].to_numpy(float)
    return x


def data(versus, col='p95', obs_ecdf='rank2d'):
    g = df[(df.role == 'data') & (df.method == 'expected_exact') & (df.versus == versus) & (df.obs_ecdf == obs_ecdf)]
    return float(g[col].iloc[0])


out = {}
rs = np.random.RandomState(0)
for label, (dv, role, ver) in {'expected_split_A': ('bigA', 'nullB', 'bigA'),
                               'expected_split_B': ('bigB', 'nullA', 'bigB'),
                               'expected_shared_A': ('bigA', 'nullA', 'bigA'),
                               'expected_shared_B': ('bigB', 'nullB', 'bigB')}.items():
    for col in ('p95', 'p99', 'max'):
        x = null(role, 'expected_exact', ver, col=col)
        d = data(dv, col)
        m, s1, s0 = x.mean(), x.std(ddof=1), x.std(ddof=0)
        # bootstrap uncertainty of the 3-sigma limit over the realisations
        lim = [np.mean(b) + 3 * np.std(b) for b in (rs.choice(x, len(x)) for _ in range(2000))]
        out[f'{label}_{col}'] = dict(data=d, null_mean=m, null_std_ddof1=s1, limit_3sig_ddof0=m + 3 * s0,
                                     limit_3sig_boot_sd=float(np.std(lim)), sep_sigma=(d - m) / s1,
                                     ratio_to_limit=d / (m + 3 * s0), null_max=float(x.max()), n=len(x),
                                     passes=bool(d <= m + 3 * s0))
for k, v in out.items():
    print(f"{k:24s} data {v['data']:.5f}  null {v['null_mean']:.5f} +- {v['null_std_ddof1']:.5f} (n={v['n']}, "
          f"max {v['null_max']:.5f})  3sig {v['limit_3sig_ddof0']:.5f} (+-{v['limit_3sig_boot_sd']:.5f})  "
          f"sep {v['sep_sigma']:.1f} sig  data/limit {v['ratio_to_limit']:.2f}  pass {v['passes']}")

# N vs N references, within this run (exact and interpolated)
print('\nN-vs-N split null, nullB (half b) realisations vs half-a references:')
nvn = []
for ref in ('bigA_m0', 'bigA_m1', 'bigA_m2', 'bigA_m3'):
    x = null('nullB', 'NvN_exact', ref=ref)
    nvn.append(x)
    print(f'  exact  ref {ref}: {x.mean():.5f} +- {x.std(ddof=1):.5f}')
for ref in ('bigA_m0', 'bigA_m1'):
    x = null('nullB', 'NvN_interp', ref=ref)
    print(f'  interp ref {ref}: {x.mean():.5f} +- {x.std(ddof=1):.5f}')
allx = np.concatenate(nvn)
print(f'  exact, pooled over 4 refs: {allx.mean():.5f} +- {allx.std(ddof=1):.5f}; '
      f'ref-to-ref sd of the mean {np.std([x.mean() for x in nvn], ddof=1):.5f}')
dm = df[(df.role == 'data') & (df.method == 'NvN_exact') & (df.versus == 'bigA_member')].p95
di = df[(df.role == 'data') & (df.method == 'NvN_interp') & (df.versus == 'bigA_member')].p95
print(f'  data vs half-a N-star members: exact {dm.mean():.5f} +- {dm.std(ddof=1):.5f}, interp {di.mean():.5f} +- {di.std(ddof=1):.5f}')
print(f'  sep (exact, pooled refs): {(dm.mean() - allx.mean()) / allx.std(ddof=1):.1f} sig')

# variance decomposition (quadrature, approximate)
sN = null('nullA', 'expected_exact', 'bigA').mean()
sN2 = null('nullB', 'expected_exact', 'bigB').mean()
sAB = null('nullB', 'expected_exact', 'bigA').mean()
sBA = null('nullA', 'expected_exact', 'bigB').mean()
sF = null('nullB', 'expected_exact', 'bigF').mean()
samp = np.mean([sN, sN2])
split = np.mean([sAB, sBA])
tex_split = np.sqrt(split ** 2 - samp ** 2)
tex_full = np.sqrt(max(sF ** 2 - samp ** 2, 0))
print(f'\nsampling only (same AST half): {sN:.5f}, {sN2:.5f} -> {samp:.5f}')
print(f'split halves (expected): {sAB:.5f}, {sBA:.5f} -> {split:.5f}; AST-half term sqrt(split^2 - samp^2) = {tex_split:.5f}')
print(f'half-b realisation vs full-list model: {sF:.5f}; implied term {tex_full:.5f} (half of the split term would be {tex_split / 2:.5f})')
print(f'quadrature check, data vs N-star (exact) {dm.mean():.5f} -> minus one sampling term {np.sqrt(dm.mean() ** 2 - samp ** 2):.5f} '
      f'(measured vs expected model {data("bigA"):.5f})')
print(f'quadrature check, N-vs-N split null predicted sqrt(2 samp^2 + tex^2) = {np.sqrt(2 * samp ** 2 + tex_split ** 2):.5f} '
      f'(measured {allx.mean():.5f}; established 0.0089)')
out['decomposition'] = dict(sampling=samp, split_expected=split, ast_half_term=tex_split, half_b_vs_full=sF,
                            implied_full_list_term=tex_full)
json.dump(out, open(f'{OUT}/summary.json', 'w'), indent=1)
