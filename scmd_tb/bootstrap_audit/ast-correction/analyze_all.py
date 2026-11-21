import glob, numpy as np, pandas as pd
OUT = '/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/ast-correction'
G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
BEST, DATA6 = 0.0168, 0.0178
fs = sorted(glob.glob(OUT + '/results_*.csv') + glob.glob(OUT + '/results2_*.csv') + glob.glob(OUT + '/results3_*.csv'))
df = pd.concat([pd.read_csv(f, dtype={'tag': str}) for f in fs], ignore_index=True)
df = df[~df.design.str.startswith('verify')].copy()
# reference identity: ('F'|'A'|'S'|'SA', split, ref seed); table refs 'T'
def refid(r):
    t = r.tag
    if t == 'T':
        return 'T-full' if r.design == 'shared' else 'T-A0'
    if 'r' in t:                       # '0r1000' -> split 0, ref seed 1000
        s, q = t.split('r')
    else:
        s, q = t, str(1000 + int(t))
    kind = {'shared': 'F', 'rmirror': 'F', 'astboot': 'F', 'mirror': f'A{s}', 'split': f'A{s}',
            'smoothed': 'S', 'smoothedA': f'SA{s}', 'split_fixedcentre': f'Afix{s}'}[r.design]
    return f'{kind}_{q}'
df['ref'] = df.apply(refid, axis=1)
# realisation list identity
def listid(r):
    t = r.tag
    s = '0' if t == 'T' else (t.split('r')[0] if 'r' in t else t)
    return {'shared': 'full', 'mirror': 'full', 'smoothed': 'full', 'smoothedA': 'full',
            'split': f'B{s}', 'rmirror': f'B{s}', 'astboot': f'L*{s}', 'split_fixedcentre': f'Bfix{s}'}[r.design]
df['list'] = df.apply(listid, axis=1)
n0 = len(df)
df = df.drop_duplicates(['design', 'ref', 'list', 'seed'])      # shared F_1005 appears as tag '5' and '0r1005'
print(f'{n0} rows, {len(df)} after removing duplicate (design, ref, list, seed)')
df.to_csv(OUT + '/all_results.csv', index=False)

print('\n=== per (design, reference, realisation list)')
g = df.groupby(['design', 'ref', 'list']).ks
t = pd.DataFrame({'n': g.size(), 'mean': g.mean(), 'std': g.std(ddof=0), 'min': g.min(), 'max': g.max()})
print(t.round(5).to_string())

print('\n=== marginal (pooled) nulls: every realisation of every reference / AST split counts once')
rows = []
for d, x in df.groupby('design'):
    gm = x.groupby(['ref', 'list']).ks
    m, s = x.ks.mean(), x.ks.std(ddof=0)
    rows.append(dict(design=d, n=len(x), combos=gm.ngroups, mean=m, std=s, within=np.sqrt(gm.var(ddof=0).mean()),
                     between=gm.mean().std(ddof=1), sig3=m + 3 * s, q9973=np.percentile(x.ks, 99.73), max=x.ks.max(),
                     z_best=(BEST - m) / s))
P = pd.DataFrame(rows).set_index('design')
print(P.round(5).to_string())
mir = df[df.design.isin(['mirror', 'rmirror', 'astboot'])]
m, s = mir.ks.mean(), mir.ks.std(ddof=0)
print(f'\nmirror-type pooled (mirror + rmirror + astboot): n={len(mir)} combos={mir.groupby(["design", "ref", "list"]).ngroups} '
      f'mean {m:.5f} std {s:.5f} 3sig {m + 3 * s:.5f} q99.73 {np.percentile(mir.ks, 99.73):.5f} max {mir.ks.max():.5f}; '
      f'best 0.0168 at z={(BEST - m) / s:.1f}; data-vs-6-reals 0.0178 at z={(DATA6 - m) / s:.1f}')
cm = mir.groupby(['design', 'ref', 'list']).ks.mean()
print(f'   combo means: min {cm.min():.5f} median {cm.median():.5f} max {cm.max():.5f}; between-combo sd {cm.std(ddof=1):.5f}')
mt = df[(df.design == 'mirror') & (df.ref == 'T-A0')].ks
print(f'mirror with the table split + table reference (T-A0 vs full list): n={len(mt)} mean {mt.mean():.5f} std {mt.std(ddof=0):.5f} '
      f'3sig {mt.mean() + 3 * mt.std(ddof=0):.5f}; best at z={(BEST - mt.mean()) / mt.std(ddof=0):.1f}')

print('\n=== split design: the table split (seed 0) vs other splits')
sp = df[df.design == 'split'].groupby(['ref', 'list']).ks.agg(['size', 'mean', 'std'])
print(sp.round(5).to_string())
fc = df[df.design == 'split_fixedcentre'].groupby(['ref', 'list']).ks.agg(['size', 'mean', 'std'])
print('fixed-centre halves:'); print(fc.round(5).to_string())
for q in ('1000', '1005'):
    a = df[(df.design == 'split') & (df.ref == f'A0_{q}')].set_index('seed').ks
    b = df[(df.design == 'split_fixedcentre') & (df.ref == f'Afix0_{q}')].set_index('seed').ks
    j = a.index.intersection(b.index)
    print(f'  split seed 0, ref {q}: per-half centres {a.loc[j].mean():.5f} vs full-list centre {b.loc[j].mean():.5f} (paired diff {(b.loc[j] - a.loc[j]).mean():+.5f}, n={len(j)})')

print('\n=== frozen-reference effect: shared design (same full-list realisations 0..K-1) by reference')
shr = df[df.design == 'shared'].groupby('ref').ks.agg(['size', 'mean', 'std'])
print(shr.round(5).T.to_string())
c = df[(df.design == 'shared') & (df.seed < 12)].pivot_table(index='seed', columns='ref', values='ks')
print(f'reference-to-reference sd of the shared-null mean (seeds 0..11, {c.shape[1]} refs): {c.mean().std(ddof=1):.5f}; '
      f'mean within-ref sd {c.std(ddof=0).mean():.5f}')

print('\n=== lambda test pooled over matched (split, ref seed) combos: (E[KS^2]mirror - E[KS^2]shared)/(E[KS^2]split - E[KS^2]shared), 1/4 expected')
def e2(d, ref, lst):
    x = df[(df.design == d) & (df.ref == ref) & (df.list == lst)].ks
    return (x ** 2).mean() if len(x) else np.nan
num = den = numr = 0.0
for ref in sorted(df[df.design == 'mirror'].ref.unique()):
    if ref == 'T-A0':
        es, em, ep, er = e2('shared', 'T-full', 'full'), e2('mirror', ref, 'full'), e2('split', ref, 'B0'), np.nan
    else:
        s, q = ref[1:].split('_')
        es, em, ep, er = e2('shared', f'F_{q}', 'full'), e2('mirror', ref, 'full'), e2('split', ref, f'B{s}'), e2('rmirror', f'F_{q}', f'B{s}')
    print(f'  {ref:8s}: shared {es:.2e} mirror {em:.2e} split {ep:.2e} rmirror {er:.2e} -> mirror ratio {(em - es) / (ep - es):5.2f}  rmirror ratio {(er - es) / (ep - es):5.2f}')
    num += em - es; den += ep - es
    if np.isfinite(er): numr += er - es
print(f'  pooled mirror ratio {num / den:.2f}')

print('\n=== smoothed error model (task suggestion) vs controls, paired by realisation')
for a, b in [('smoothed', 'shared'), ('smoothedA', 'mirror'), ('smoothedA', 'smoothed'), ('mirror', 'shared')]:
    x = df[df.design == a]; y = df[df.design == b]
    x = x.assign(k=x.ref.str.split('_').str[-1]).set_index(['k', 'seed']).ks
    y = y.assign(k=y.ref.str.split('_').str[-1]).set_index(['k', 'seed']).ks
    y = y[~y.index.duplicated()]
    j = x.index.intersection(y.index)
    j = [i for i in j if i[0] in ('1001', '1002')]
    dd = x.loc[j].values - y.loc[j].values
    print(f'  {a} - {b} (refs 1001, 1002): n={len(dd)} mean {dd.mean():+.5f} SE {dd.std(ddof=1) / np.sqrt(len(dd)):.5f}')
