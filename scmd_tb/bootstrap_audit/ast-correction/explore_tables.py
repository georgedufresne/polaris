import numpy as np, pandas as pd
G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
sh = pd.read_csv(G + '/bootstrap_2d_ks_results.csv')
sp = pd.read_csv(G + '_split/bootstrap_2d_ks_results.csv')
for name, t in [('shared', sh), ('split', sp)]:
    k = t.ks_statistic.values
    print(name, len(k), 'mean %.6f std %.6f median %.6f  3sig %.6f  min %.6f max %.6f' % (k.mean(), k.std(), np.median(k), k.mean() + 3 * k.std(), k.min(), k.max()),
          'n_obs', t.n_observed.unique(), 'n_sim', np.unique(t.n_simulated)[:5], 'ids', t.bootstrap_number.min(), t.bootstrap_number.max(), t.bootstrap_number.nunique())
    print('   mean(KS^2) %.4e  rms %.6f  skew %.3f' % ((k**2).mean(), np.sqrt((k**2).mean()), pd.Series(k).skew()))
# paired by bootstrap number (same realisation seed, different AST list for realisations though)
m = sh.merge(sp, on='bootstrap_number', suffixes=('_sh', '_sp'))
print('paired rows', len(m), 'corr(KS_shared, KS_split) = %.3f' % np.corrcoef(m.ks_statistic_sh, m.ks_statistic_sp)[0, 1])
print('spearman %.3f' % m[['ks_statistic_sh', 'ks_statistic_sp']].corr('spearman').iloc[0, 1])
