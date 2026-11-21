#!/usr/bin/env python3
"""
Uncertainty budget for the sCMD-only posteriors, OG vs NEW.

Three separate things get called "the uncertainty" and they are not the same:

  1. POSTERIOR WIDTH (sigma, and the asymmetric P16/P84 interval) -- the
     uncertainty ON THE PARAMETER.  This is what you quote as "8.53 +/- 0.64 Gyr".
     It is what analyze_scmd_only.py already reports as 'Std'.

  2. SAMPLING ERROR ON THE POSTERIOR MEAN -- how well the finite grid of MC sets
     pins down the centre of that posterior.  Rows are NOT independent: all the
     (age, DM, EBV) rows from one MC set share the same physics, so the
     independent unit is the MC SET, not the row.  Reported here as an effective
     sample size clustered on simulation_id, plus a cluster bootstrap.

  3. UNCERTAINTY ON THE OG->NEW SHIFT.  The two grids reuse the SAME var-file
     realisations (same file_id -> same FeH/He/alphafe/CMIXLA), so the two
     samples are PAIRED, not independent.  The paired difference is far better
     determined than sqrt(SE_og^2 + SE_new^2) would suggest.  Computed on the MC
     sets common to both, resampling sets (not rows) with replacement.

NOTE: this is a cluster bootstrap over rows that already exist -- it is NOT a
re-run of the KS null-distribution bootstrap, which the user declined and which
would require recomputing sCMDs.  Cost here is a few seconds.
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

OG_HOME  = os.environ.get('OG_HOME',  '/dartfs-hpc/rc/home/8/f007058')
NEW_HOME = os.environ.get('NEW_HOME', '/dartfs-hpc/scratch/f007058/analysis_noco')
BOOTSTRAP = f'{OG_HOME}/Bootstrap/bootstrap_2d_ks_results.csv'
INITIAL_VARS = f'{OG_HOME}/initial_vars.csv'
OUT = os.environ.get('UNC_OUT', '/dartfs-hpc/scratch/f007058/scmd_og_vs_noco_uncertainties.txt')

N_BOOT = 2000
RNG = np.random.default_rng(20261010)

# --- binning, copied verbatim from analyze_scmd_only.py -----------------------
age_bins  = np.arange(3.95, 12.15, 0.1)
feh_bins  = np.arange(-0.505, 1.005, 0.02)
dm_bins   = np.arange(10.005, 16.005, 0.02)
ebv_bins  = np.arange(-0.005, 0.805, 0.01)
he_bins   = np.arange(0.195, 0.405, 0.002)
afe_bins  = np.array([-0.3, -0.1, 0.1, 0.3, 0.5])
cmix_bins = np.arange(0.995, 2.505, 0.01)
ctr = lambda b: (b[:-1] + b[1:]) / 2

PARAMS = [
    ('Age (Gyr)',     'age_gyr',           age_bins,  False),
    ('[Fe/H]',        'FeH',               feh_bins,  False),
    ('Dist. Mod.',    'distance_modulus',  dm_bins,   False),
    ('E(B-V)',        'extinction_ebv',    ebv_bins,  False),
    ('He (Y)',        'He Abundance',      he_bins,   False),
    ('[alpha/Fe]',    'alphafe',           afe_bins,  True),
    ('Mixing Length', 'CMIXLA',            cmix_bins, False),
]
AFE_CENTERS = np.array([-0.2, 0.0, 0.2, 0.4])


def compute_stats(vals, prob, discrete=False):
    """Identical to analyze_scmd_only.py so the numbers line up exactly."""
    prob = prob / np.sum(prob)
    mean = np.sum(vals * prob)
    cum = np.cumsum(prob)

    def interp_pct(pct):
        idx = np.searchsorted(cum, pct)
        if idx == 0:
            return vals[0]
        if idx >= len(vals):
            return vals[-1]
        if discrete:
            return vals[idx]
        frac = (pct - cum[idx - 1]) / (cum[idx] - cum[idx - 1] + 1e-10)
        return vals[idx - 1] + frac * (vals[idx] - vals[idx - 1])

    return {
        'mean': mean,
        'median': interp_pct(0.50),
        'mode': vals[np.argmax(prob)],
        'p16': interp_pct(0.16),
        'p84': interp_pct(0.84),
        'std': np.sqrt(np.sum(prob * (vals - mean) ** 2)),
    }


# --- KS -> weight, shared by both sides --------------------------------------
boot = pd.read_csv(BOOTSTRAP)
ks_kde = gaussian_kde(boot['ks_statistic'].values, bw_method='scott')
initial_vars = pd.read_csv(INITIAL_VARS)
initial_vars['file_id'] = initial_vars['file_id'].astype(int)


def load(home):
    """Reproduce analyze_scmd_only.py's sample exactly, unnormalised weights."""
    d = pd.read_csv(f'{home}/filtered_ks_scmd.csv')
    d['w'] = ks_kde(d['ks_statistic'].values)
    # keep only the best (lowest KS) DM+EBV combo per (MC set, age)
    d = d.loc[d.groupby(['simulation_id', 'age_gyr'])['ks_statistic'].idxmin()]
    d = pd.merge(d, initial_vars, left_on='simulation_id', right_on='file_id', how='left')
    return d.reset_index(drop=True)


def posterior_mean(vals, w, bins, discrete):
    centers = AFE_CENTERS if discrete else ctr(bins)
    h, _ = np.histogram(vals, bins=bins, weights=w)
    if h.sum() <= 0:
        return np.nan
    return compute_stats(centers, h / h.sum(), discrete)['mean']


def full_stats(d, col, bins, discrete):
    m = d[col].notna().values
    centers = AFE_CENTERS if discrete else ctr(bins)
    h, _ = np.histogram(d.loc[m, col].values, bins=bins, weights=d.loc[m, 'w'].values)
    return compute_stats(centers, h / h.sum(), discrete)


def set_index_map(d):
    """row indices grouped by MC set, for cluster resampling."""
    g = d.groupby('simulation_id').indices
    return list(g.keys()), g


def cluster_boot_mean(d, col, bins, discrete, n=N_BOOT, sets=None, groups=None):
    """Bootstrap the posterior mean, resampling MC SETS with replacement."""
    if sets is None:
        sets, groups = set_index_map(d)
    sets = np.asarray(sets)
    vals_all = d[col].values
    w_all = d['w'].values
    out = np.empty(n)
    for i in range(n):
        pick = RNG.integers(0, len(sets), len(sets))
        idx = np.concatenate([groups[sets[j]] for j in pick])
        v, w = vals_all[idx], w_all[idx]
        ok = ~np.isnan(v)
        out[i] = posterior_mean(v[ok], w[ok], bins, discrete)
    return out


def ess_clustered(d):
    """Effective sample size with the MC set as the independent unit."""
    W = d.groupby('simulation_id')['w'].sum().values
    return W.sum() ** 2 / np.sum(W ** 2)


og, new = load(OG_HOME), load(NEW_HOME)
common = np.intersect1d(og['simulation_id'].unique(), new['simulation_id'].unique())
og_c = og[og['simulation_id'].isin(common)].reset_index(drop=True)
new_c = new[new['simulation_id'].isin(common)].reset_index(drop=True)

ess_og, ess_new = ess_clustered(og), ess_clustered(new)

L = []
A = L.append
A('=' * 100)
A('UNCERTAINTY BUDGET FOR THE sCMD-ONLY POSTERIORS  --  ORIGINAL vs NEW (fixed physics, no overshoot)')
A('=' * 100)
A('')
A(f'  rows (best DM+EBV per set+age) :  OG {len(og):>6,}      NEW {len(new):>6,}')
A(f'  distinct MC sets               :  OG {og["simulation_id"].nunique():>6,}      NEW {new["simulation_id"].nunique():>6,}')
A(f'  MC sets common to both         :  {len(common):,}  <- the paired sample')
A('')
A('  Rows are not independent: every (age, DM, EBV) row inside one MC set shares the same')
A('  physics realisation, so the independent unit is the MC SET.  Effective sample size,')
A('  clustered on simulation_id and weighted by P_KS:')
A(f'    OG  : ESS = {ess_og:7.1f} independent MC sets (of {og["simulation_id"].nunique()} present)')
A(f'    NEW : ESS = {ess_new:7.1f} independent MC sets (of {new["simulation_id"].nunique()} present)')
A('  The KS weights are concentrated, so the grid behaves like far fewer sets than it holds.')
A('')

rows = []
A('-' * 100)
A('1. UNCERTAINTY ON THE PARAMETER  (posterior width -- this is what you quote)')
A('-' * 100)
A(f'{"Parameter":<16}{"ORIGINAL":>34}{"NEW":>34}')
A(f'{"":<16}{"mean +/- sigma":>22}{"68% asym":>12}{"mean +/- sigma":>22}{"68% asym":>12}')
for name, col, bins, disc in PARAMS:
    so, sn = full_stats(og, col, bins, disc), full_stats(new, col, bins, disc)
    rows.append((name, col, bins, disc, so, sn))
    f = 4 if name in ('He (Y)',) else 3
    qo = f'{so["mean"]:.{f}f} +/- {so["std"]:.{f}f}'
    qn = f'{sn["mean"]:.{f}f} +/- {sn["std"]:.{f}f}'
    ao = f'+{so["p84"]-so["median"]:.{f}f}/-{so["median"]-so["p16"]:.{f}f}'
    an = f'+{sn["p84"]-sn["median"]:.{f}f}/-{sn["median"]-sn["p16"]:.{f}f}'
    A(f'{name:<16}{qo:>22}{ao:>12}{qn:>22}{an:>12}')
A('')
A('  "68% asym" is P16/P84 about the MEDIAN -- use it when the posterior is skewed,')
A('  which [alpha/Fe] badly is (median 0.000, P84 0.200).')
A('')

A('-' * 100)
A('2. SAMPLING ERROR ON THE POSTERIOR MEAN  (how well the grid pins down the centre)')
A('-' * 100)
A(f'{"Parameter":<16}{"OG sigma":>12}{"OG SE":>12}{"NEW sigma":>12}{"NEW SE":>12}{"SE as % of sigma":>20}')
se = {}
for name, col, bins, disc, so, sn in rows:
    bo = cluster_boot_mean(og, col, bins, disc)
    bn = cluster_boot_mean(new, col, bins, disc)
    se[name] = (np.std(bo, ddof=1), np.std(bn, ddof=1))
    f = 4 if name == 'He (Y)' else 3
    pct = f'{se[name][0]/so["std"]*100:.0f}% / {se[name][1]/sn["std"]*100:.0f}%'
    A(f'{name:<16}{so["std"]:>12.{f}f}{se[name][0]:>12.{f}f}{sn["std"]:>12.{f}f}{se[name][1]:>12.{f}f}{pct:>20}')
A('')
A(f'  Cluster bootstrap, {N_BOOT} resamples of MC sets with replacement.')
A('')

A('-' * 100)
A('3. IS THE OG -> NEW SHIFT SIGNIFICANT?  (paired, on the MC sets common to both)')
A('-' * 100)
A(f'{"Parameter":<16}{"shift":>12}{"SE(shift)":>12}{"shift/SE":>12}{"95% CI on the shift":>30}')
sig = []
sets_c = np.asarray(sorted(common))
g_og = og_c.groupby('simulation_id').indices
g_new = new_c.groupby('simulation_id').indices
for name, col, bins, disc, _so, _sn in rows:
    vo, wo = og_c[col].values, og_c['w'].values
    vn, wn = new_c[col].values, new_c['w'].values
    m_og = posterior_mean(*(lambda v, w: (v[~np.isnan(v)], w[~np.isnan(v)]))(vo, wo), bins, disc)
    m_new = posterior_mean(*(lambda v, w: (v[~np.isnan(v)], w[~np.isnan(v)]))(vn, wn), bins, disc)
    shift = m_new - m_og
    diffs = np.empty(N_BOOT)
    for i in range(N_BOOT):
        pick = RNG.integers(0, len(sets_c), len(sets_c))
        io = np.concatenate([g_og[sets_c[j]] for j in pick])
        inw = np.concatenate([g_new[sets_c[j]] for j in pick])
        a, b = vo[io], wo[io]
        ok = ~np.isnan(a)
        ma = posterior_mean(a[ok], b[ok], bins, disc)
        a2, b2 = vn[inw], wn[inw]
        ok2 = ~np.isnan(a2)
        diffs[i] = posterior_mean(a2[ok2], b2[ok2], bins, disc) - ma
    s = np.std(diffs, ddof=1)
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    z = shift / s if s > 0 else np.nan
    sig.append((name, shift, s, z, lo, hi))
    f = 4 if name == 'He (Y)' else 3
    ci = f'[{lo:+.{f}f}, {hi:+.{f}f}]'
    A(f'{name:<16}{shift:>+12.{f}f}{s:>12.{f}f}{z:>12.2f}{ci:>30}')
A('')
A('  shift/SE is a z-score: |z| < 2 means the change is consistent with zero.')
A('')
A('-' * 100)
A('VERDICT')
A('-' * 100)
for name, shift, s, z, lo, hi in sig:
    tag = 'SIGNIFICANT' if abs(z) >= 2 else 'consistent with zero'
    A(f'  {name:<16} z = {z:+5.2f}   {tag}')
A('')
A('=' * 100)

txt = '\n'.join(L)
open(OUT, 'w').write(txt + '\n')
print(txt)
print(f'\nSaved -> {OUT}')
