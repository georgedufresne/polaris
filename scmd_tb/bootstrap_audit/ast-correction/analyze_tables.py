import numpy as np, pandas as pd
G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
sh = pd.read_csv(G + '/bootstrap_2d_ks_results.csv').sort_values('bootstrap_number').ks_statistic.values
sp = pd.read_csv(G + '_split/bootstrap_2d_ks_results.csv').sort_values('bootstrap_number').ks_statistic.values
BEST = 0.0168
def desc(name, k):
    m, s = k.mean(), k.std()
    print(f'{name:38s} mean {m:.5f} std {s:.5f} 3sig {m + 3 * s:.5f} | median {np.median(k):.5f} q99.73 {np.percentile(k, 99.73):.5f} max {k.max():.5f} | best z = {(BEST - m) / s:5.1f}')
    return m, s
print('--- tables')
desc('shared (10k)', sh); desc('split (10k)', sp)
print('--- mirror estimates, null^2 = shared^2 + (split^2 - shared^2)/4 = 0.75 shared^2 + 0.25 split^2')
m1 = np.sqrt(0.75 * sh.mean() ** 2 + 0.25 * sp.mean() ** 2)
m2 = np.sqrt(0.75 * (sh ** 2).mean() + 0.25 * (sp ** 2).mean())
m3 = np.sqrt(0.75 * np.median(sh) ** 2 + 0.25 * np.median(sp) ** 2)
print(f'on means {m1:.5f}; on mean squares (rms) {m2:.5f}; on medians {m3:.5f}')
qm = np.sqrt(0.75 * np.sort(sh) ** 2 + 0.25 * np.sort(sp) ** 2)      # quantile (rank) matched
mq = desc('quantile-matched quadrature', qm)
pm = np.sqrt(0.75 * sh ** 2 + 0.25 * sp ** 2)                        # paired by realisation seed
desc('paired-by-seed quadrature', pm)
ql = np.sort(sh) + 0.5 * (np.sort(sp) - np.sort(sh))                  # linear alternative (lambda=1/2 in KS)
desc('quantile-matched linear midpoint', ql)
print('--- additivity check: if the AST term added a constant in quadrature, split_q^2 - shared_q^2 would be flat in p')
ps = [0.01, 0.05, 0.16, 0.25, 0.5, 0.75, 0.84, 0.95, 0.99]
qs, qp = np.quantile(sh, ps), np.quantile(sp, ps)
print(' p     shared_q  split_q   split_q^2-shared_q^2   split_q-shared_q   ratio')
for p, a, b in zip(ps, qs, qp):
    print(f'{p:4.2f}  {a:.5f}  {b:.5f}   {b * b - a * a:.3e}            {b - a:.5f}          {b / a:.2f}')
print('constant-quadrature model: std(split) would be', end=' ')
A2 = (sp ** 2).mean() - (sh ** 2).mean()
pred = np.sqrt(np.sort(sh) ** 2 + A2)
print(f'{pred.std():.5f} (observed {sp.std():.5f}); predicted split quantiles 5/50/95% {np.quantile(pred, [.05, .5, .95]).round(5)} vs observed {np.quantile(sp, [.05, .5, .95]).round(5)}')
print('corr(shared_i, split_i) paired by seed = %.3f' % np.corrcoef(sh, sp)[0, 1])
# standard error of the table means and of the formula, from the per-row values
n = len(sh)
print(f'SE of table means: shared {sh.std() / np.sqrt(n):.6f}, split {sp.std() / np.sqrt(n):.6f}')
