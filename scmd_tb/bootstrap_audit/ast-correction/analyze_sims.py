import glob, numpy as np, pandas as pd
OUT = '/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/ast-correction'
G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
BEST = 0.0168
df = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(OUT + '/results_*.csv'))], ignore_index=True)
df['tag'] = df['tag'].astype(str)
print('rows', len(df), 'n_sim range', df.n_sim.min(), df.n_sim.max(), 'n_ref range', df.n_ref.min(), df.n_ref.max())
sh = pd.read_csv(G + '/bootstrap_2d_ks_results.csv').set_index('bootstrap_number')['ks_statistic']
sp = pd.read_csv(G + '_split/bootstrap_2d_ks_results.csv').set_index('bootstrap_number')['ks_statistic']

v = df[df.role == 'V']
if len(v):
    print('\n=== V: reproduction of the 10,000-row tables (same references, same seeds)')
    for d, tab in [('verify_split', sp), ('verify_shared', sh), ('split', sp), ('shared', sh)]:
        x = v[(v.design == d)]
        if len(x):
            diff = x.ks.values - tab.loc[x.seed.values].values
            print(f'  {d:14s} n={len(x):2d}  max |code - table| = {np.abs(diff).max():.2e}')

print('\n=== per design and tag (tag = split seed s / resample b; T = table references)')
g = df[~df.design.str.startswith('verify')].groupby(['design', 'tag']).ks
t = g.agg(['count', 'mean', 'std', 'median', 'min', 'max'])
t['std'] = g.std(ddof=0)
print(t.round(5).to_string())

print('\n=== pooled over tags (fresh refs s=1..4 only; T listed separately)')
fresh = df[(df.tag != 'T') & ~df.design.str.startswith('verify')]
rows = []
for d, x in fresh.groupby('design'):
    m_tag = x.groupby('tag').ks.mean()
    rows.append(dict(design=d, n=len(x), ntag=x.tag.nunique(), mean=x.ks.mean(), std_pooled=x.ks.std(ddof=0),
                     within_std=np.sqrt(x.groupby('tag').ks.var(ddof=0).mean()), between_std_of_tag_means=m_tag.std(ddof=1),
                     rms=np.sqrt((x.ks ** 2).mean())))
P = pd.DataFrame(rows).set_index('design')
P['3sig_pooled'] = P['mean'] + 3 * P['std_pooled']
P['best_z_pooled'] = (BEST - P['mean']) / P['std_pooled']
print(P.round(5).to_string())

print('\n=== paired differences (same realisation sCMD, references with the same seed)')
def paired(a, b, la, lb):
    x = fresh[fresh.design == a].set_index(['tag', 'seed']).ks
    y = fresh[fresh.design == b].set_index(['tag', 'seed']).ks
    j = x.index.intersection(y.index)
    dd = x.loc[j] - y.loc[j]
    print(f'  {la} - {lb}: n={len(j)} mean {dd.mean():+.5f} (SE {dd.std(ddof=1) / np.sqrt(len(j)):.5f}); by tag:',
          ' '.join(f'{s}:{dd.xs(s, level=0).mean():+.5f}' for s in sorted(set(j.get_level_values(0)))))
paired('mirror', 'shared', 'mirror(A ref, full reals)', 'shared(F ref, full reals)')
paired('split', 'rmirror', 'split(A ref, B reals)', 'rmirror(F ref, B reals)')
paired('smoothed', 'shared', 'smoothed(F-smoothed ref)', 'shared')
paired('smoothedA', 'mirror', 'smoothedA(A-smoothed ref)', 'mirror(A ref)')
paired('smoothedA', 'smoothed', 'smoothedA', 'smoothed')
paired('astboot', 'shared', 'astboot(F ref, resampled-list reals)', 'shared(F ref, full reals)')

print('\n=== lambda test per tag: is (E[KS^2]mirror - E[KS^2]shared) = (E[KS^2]split - E[KS^2]shared)/4 ?')
for s in sorted(set(fresh.tag)):
    e = {d: (fresh[(fresh.design == d) & (fresh.tag == s)].ks ** 2).mean() for d in ('shared', 'mirror', 'split', 'rmirror')}
    if all(np.isfinite(list(e.values()))):
        dsp = e['split'] - e['shared']
        print(f'  s={s}: E[KS^2] shared {e["shared"]:.2e} mirror {e["mirror"]:.2e} rmirror {e["rmirror"]:.2e} split {e["split"]:.2e} | '
              f'(mirror-shared)/(split-shared) = {(e["mirror"] - e["shared"]) / dsp:.2f}, (rmirror-shared)/(split-shared) = {(e["rmirror"] - e["shared"]) / dsp:.2f}  (1/4 expected)')
T = df[(df.tag == 'T') & ~df.design.str.startswith('verify')]
if len(T):
    e = {d: (T[T.design == d].ks ** 2).mean() for d in ('shared', 'mirror', 'split')}
    print(f'  T (table refs): E[KS^2] shared {e["shared"]:.2e} mirror {e["mirror"]:.2e} split {e["split"]:.2e} | ratio {(e["mirror"] - e["shared"]) / (e["split"] - e["shared"]):.2f}')
    for d in ('shared', 'mirror', 'split'):
        x = T[T.design == d].ks
        print(f'    T {d:7s}: n={len(x)} mean {x.mean():.5f} std {x.std(ddof=0):.5f}')
df.to_csv(OUT + '/all_results.csv', index=False)
