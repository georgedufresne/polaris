import glob, numpy as np, pandas as pd
O = '/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/reference-selection/out'
sp = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(O + '/phase2_split_*.csv'))]).sort_values('seed').reset_index(drop=True)
sh = pd.concat([pd.read_csv(f) for f in sorted(glob.glob(O + '/phase2_shared_*.csv'))]).sort_values('seed').reset_index(drop=True)
pub_split = pd.read_csv('/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02_split/bootstrap_2d_ks_results.csv').set_index('bootstrap_number')
pub_shared = pd.read_csv('/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02/bootstrap_2d_ks_results.csv').set_index('bootstrap_number')
print('n realisations split', len(sp), 'shared', len(sh))
d1 = sp.R0_ks.values - pub_split.loc[sp.seed, 'ks_statistic'].values
d2 = sh.S0_ks.values - pub_shared.loc[sh.seed, 'ks_statistic'].values
print(f'VALIDATION R0 vs published split null (same seeds): max |diff| {np.abs(d1).max():.2e}, mean diff {d1.mean():.2e}')
print(f'VALIDATION S0 vs published shared null (same seeds): max |diff| {np.abs(d2).max():.2e}, mean diff {d2.mean():.2e}')
print(f'published split null seeds 0-47: mean {pub_split.loc[sp.seed].ks_statistic.mean():.6f}; all 10000: mean {pub_split.ks_statistic.mean():.6f} std {pub_split.ks_statistic.std(ddof=0):.6f}')
print(f'published shared null seeds 0-47: mean {pub_shared.loc[sh.seed].ks_statistic.mean():.6f}; all 10000: mean {pub_shared.ks_statistic.mean():.6f} std {pub_shared.ks_statistic.std(ddof=0):.6f}')

def twoway(Y, label):
    # Y: J realisations x K references; random-effects two-way ANOVA without replication
    J, K = Y.shape
    gm = Y.mean()
    rm = Y.mean(axis=0)  # per reference
    cm = Y.mean(axis=1)  # per realisation
    ss_ref = J * ((rm - gm) ** 2).sum(); ss_real = K * ((cm - gm) ** 2).sum()
    ss_res = ((Y - rm[None, :] - cm[:, None] + gm) ** 2).sum()
    ms_ref = ss_ref / (K - 1); ms_real = ss_real / (J - 1); ms_res = ss_res / ((J - 1) * (K - 1))
    sB2 = max((ms_ref - ms_res) / J, 0); sM2 = max((ms_real - ms_res) / K, 0); sE2 = ms_res
    F = ms_ref / ms_res
    from scipy.stats import f as fdist
    p = fdist.sf(F, K - 1, (J - 1) * (K - 1))
    print(f'\n[{label}] {J} realisations x {K} references')
    print(f'  per-reference mean (std over realisations, ddof 1):')
    for k in range(K):
        print(f'    {names[k]:>4}: mean {rm[k]:.6f}  std {Y[:, k].std(ddof=1):.6f}  mean+3sd {rm[k] + 3 * Y[:, k].std(ddof=1):.6f}  min {Y[:, k].min():.6f} max {Y[:, k].max():.6f}')
    print(f'  spread of the {K} reference means: std {rm.std(ddof=1):.6f}, range {rm.min():.6f}-{rm.max():.6f}; expected from realisation noise alone sqrt(MS_res/J) = {np.sqrt(ms_res / J):.6f}')
    print(f'  F(reference) = {F:.1f}, p = {p:.2e}')
    print(f'  variance components: sigma_ref {np.sqrt(sB2):.6f}, sigma_realisation {np.sqrt(sM2):.6f}, sigma_interaction {np.sqrt(sE2):.6f}')
    cond = np.sqrt(sM2 + sE2); unc = np.sqrt(sB2 + sM2 + sE2)
    print(f'  conditional (one reference) sd {cond:.6f}; unconditional sd {unc:.6f}; ratio {unc / cond:.3f}')
    print(f'  pooled over all references: mean {gm:.6f}, sd {Y.std(ddof=1):.6f}, mean+3sd {gm + 3 * Y.std(ddof=1):.6f}; mean+3*unconditional sd {gm + 3 * unc:.6f}')
    # mean+3sd per reference, range
    th = rm + 3 * Y.std(axis=0, ddof=1)
    print(f'  3-sigma limit per reference: {th.min():.6f} .. {th.max():.6f} (std {th.std(ddof=1):.6f})')
    return rm, sB2, sM2, sE2

names = [f'R{k}' for k in range(10)]
Ys = sp[[f'{n}_ks' for n in names]].values
rmS, sB2s, sM2s, sE2s = twoway(Ys, 'SPLIT: refs half_a, realisations half_b; R0 = official reference')
print(f'  R0 (official) minus mean of fresh refs R1-R9: {rmS[0] - rmS[1:].mean():+.6f}; R0 rank among 10 means: {int((rmS < rmS[0]).sum()) + 1}')
# fresh refs only
names_f = names[1:]
names_bak = names; names = names_f
twoway(Ys[:, 1:], 'SPLIT, fresh references R1-R9 only')
names = names_bak
names = ['S0', 'G0', 'S1', 'S2', 'S3', 'S4', 'S5']
Yh = sh[[f'{n}_ks' for n in names]].values
rmH, _, _, _ = twoway(Yh, 'SHARED: refs and realisations full AST list; S0 = official ref, G0 = grid code exact selected sCMD')
print(f'  S0 vs G0 (same grid seed, ridge traced differently): mean diff {np.mean(sh.S0_ks - sh.G0_ks):+.2e}, max |diff| {np.abs(sh.S0_ks - sh.G0_ks).max():.2e}')
print(f'  grid-selected refs (S0, G0) mean {rmH[:2].mean():.6f} vs fresh refs S1-S5 mean {rmH[2:].mean():.6f} (sd of fresh means {rmH[2:].std(ddof=1):.6f})')
names = ['S1', 'S2', 'S3', 'S4', 'S5']
twoway(sh[[f'{n}_ks' for n in names]].values, 'SHARED, fresh references S1-S5 only')
# data vs realisations (b)
x = sh.data_ks.values
print(f'\n(b) data vs {len(x)} fresh realisations (latest code, full AST): mean {x.mean():.6f} sd {x.std(ddof=1):.6f} sem {x.std(ddof=1)/np.sqrt(len(x)):.6f} median {np.median(x):.6f} min {x.min():.6f} max {x.max():.6f}')
print(f'    fraction <= 0.016828: {np.mean(x <= 0.016828):.3f} ({int(np.sum(x <= 0.016828))}/{len(x)}); z = {(0.016828 - x.mean()) / x.std(ddof=1):+.2f}')
gc = pd.read_csv(O + '/gridcode_data_vs_fresh.csv')
y = gc.data_ks.values
print(f'    data vs {len(y)} fresh realisations with the GRID code: mean {y.mean():.6f} sd {y.std(ddof=1):.6f} sem {y.std(ddof=1)/np.sqrt(len(y)):.6f}; fraction <= 0.016828: {np.mean(y <= 0.016828):.3f}; z = {(0.016828 - y.mean()) / y.std(ddof=1):+.2f}')
# paired comparison latest vs grid code on the same seeds 0..23
m = sh.set_index('seed').loc[gc.seed, 'data_ks'].values
print(f'    same seeds 0-23, latest minus grid code: mean {np.mean(m - y):+.6f} sd {np.std(m - y, ddof=1):.6f} (t = {np.mean(m - y) / (np.std(m - y, ddof=1) / np.sqrt(len(y))):+.2f}); corr {np.corrcoef(m, y)[0, 1]:.2f}')
# the references S1..S5 are full-AST realisations too: KS(data, S_k) computed in phase 2? (not computed: only realisation->data). 
print('\n  p99 and max |diff| for data vs realisations: p99 mean %.5f, max mean %.4f' % (sh.data_p99.mean(), sh.data_max.mean()))
sp.to_csv(O + '/split_design_all.csv', index=False); sh.to_csv(O + '/shared_design_all.csv', index=False)
