import numpy as np, pandas as pd
G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
sh = pd.read_csv(G + '/bootstrap_2d_ks_results.csv').sort_values('bootstrap_number').ks_statistic.values
sp = pd.read_csv(G + '_split/bootstrap_2d_ks_results.csv').sort_values('bootstrap_number').ks_statistic.values
# KS^2 = frozen^2 + realisation^2 + 2 x cross term; cross term has zero mean and std proportional to the
# frozen AST amplitude, so for an AST amplitude lambda (0 shared, 1/2 mirror, 1 split):
#   E[KS^2](l) = E0 + l^2 (E1 - E0),   Var[KS^2](l) = V0 + l^2 (V1 - V0)   (no cov between cross and noise terms)
E0, E1, V0, V1 = (sh**2).mean(), (sp**2).mean(), (sh**2).var(), (sp**2).var()
l = 0.5
E, V = E0 + l * l * (E1 - E0), V0 + l * l * (V1 - V0)
print(f'E[KS^2]: shared {E0:.3e} split {E1:.3e} -> mirror {E:.3e};  sd[KS^2]: shared {np.sqrt(V0):.3e} split {np.sqrt(V1):.3e} -> mirror {np.sqrt(V):.3e}')
# map back to KS with a gamma distribution for KS^2 with those moments (KS^2 is positive, skewed)
k, th = E * E / V, V / E
rng = np.random.default_rng(0)
x = np.sqrt(rng.gamma(k, th, 2_000_000))
print(f'moment-interpolated mirror (gamma on KS^2): mean {x.mean():.5f} std {x.std():.5f} 3sig {x.mean() + 3 * x.std():.5f} q99.73 {np.percentile(x, 99.73):.5f}')
# same machinery reproduces the two tables?
for nm, a in [('shared', sh), ('split', sp)]:
    e, v = (a**2).mean(), (a**2).var()
    y = np.sqrt(rng.gamma(e * e / v, v / e, 2_000_000))
    print(f'  check {nm}: gamma-model mean {y.mean():.5f} std {y.std():.5f}  vs table mean {a.mean():.5f} std {a.std():.5f}')
