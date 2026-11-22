#!/usr/bin/env python3
"""Six-way sCMD comparison for NGC 6791: the five-way (../scmd_W14_pipeline/scmd_5way_compare.py) plus W14feh = W v1.4
with [Fe/H] converted to our Z convention (gen_W14_iso.py --z-factor 1.121), so W models the same Z as OG/New.
Copy changed only in the legs, the pairs, the output name and the home inputs (scratch copies in inputs/, md5-checked).
Original docstring: Five-way sCMD comparison for NGC 6791: OG vs New vs Emul(old) vs W (v1.2, October) vs W14 (W v1.4, CONF1DENCE-W
fcfac3c). Extends ../scmd_3way_compare.py (same weights: the original bootstrap null's KDE at each KS value, one best
DM/E per (MC set, age); same paired bootstrap of the age difference on common MC sets). Each leg's parameter summary is
its own analyze_scmd_only.py output (scmd_only_parameter_summary.csv). "sets kept / fitted" = MC sets surviving the
fixed KS cut (0.0468) over MC sets with any 2D-KS result (newest combined_2d_ks_results_*.csv in the leg's analysis dir;
OG: the home copy).
KTTAU: W v1.4 was not trained on KTTAU = 1, so every leg is restricted to the MC sets with KTTAU 0 or 5
(initial_vars.csv), and the parameter posteriors are computed here from those rows (KS-weighted mean and sd) instead of
being read from each leg's scmd_only_parameter_summary.csv (which covers all sets).
Output: ../scmd_5way_comparison.txt."""
import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

S = '/dartfs-hpc/scratch/f007058'
H = '/dartfs-hpc/rc/home/8/f007058'
LEGS = {'OG': f'{S}/analysis_OG', 'New': f'{S}/analysis_noco', 'Emul(old)': f'{S}/analysis_emul',
        'W': f'{S}/analysis_W', 'W14': f'{S}/analysis_W14', 'W14feh': f'{S}/analysis_W14feh'}
FILT = {'OG': f'{S}/scmd_W14feh_pipeline/inputs/filtered_ks_scmd_OG.csv', 'New': f'{S}/analysis_noco/filtered_ks_scmd.csv',
        'Emul(old)': f'{S}/analysis_emul/filtered_ks_scmd.csv', 'W': f'{S}/analysis_W/filtered_ks_scmd.csv',
        'W14': f'{S}/analysis_W14/filtered_ks_scmd.csv', 'W14feh': f'{S}/analysis_W14feh/filtered_ks_scmd.csv'}
OUT = f'{S}/scmd_6way_comparison.txt'
N_BOOT = 2000
RNG = np.random.default_rng(20261118)

IV = pd.read_csv(f'{S}/analysis_W/initial_vars.csv').rename(columns={'file_id': 'simulation_id'}).set_index('simulation_id')
KEEP = set(IV.index[IV.KTTAU.isin([0, 5])])
boot = pd.read_csv(f'{S}/scmd_W14feh_pipeline/inputs/bootstrap_2d_ks_results.csv')
kde = gaussian_kde(boot.ks_statistic.values, bw_method='scott')
age_bins = np.arange(3.95, 12.15, 0.1)
age_ctr = (age_bins[:-1] + age_bins[1:]) / 2


def load(p):
    d = pd.read_csv(p)
    d = d[d.simulation_id.isin(KEEP)]
    d['w'] = kde(d.ks_statistic.values)
    return d.loc[d.groupby(['simulation_id', 'age_gyr'])['ks_statistic'].idxmin()].reset_index(drop=True)


def wmean(v, w, bins, ctr):
    h, _ = np.histogram(v, bins=bins, weights=w)
    if h.sum() <= 0:
        return np.nan
    p = h / h.sum()
    return float(np.sum(ctr * p))


def n_fitted(leg):
    d = LEGS[leg]
    c = sorted(glob.glob(f'{d}/combined_2d_ks_results_*.csv'), key=os.path.getmtime)
    if leg == 'OG' and not c:
        c = sorted(glob.glob(f'{H}/combined_2d_ks_results_*.csv'), key=os.path.getmtime)
    if not c:
        return np.nan
    ids = pd.read_csv(c[-1], usecols=['simulation_id']).simulation_id
    return ids[ids.isin(KEEP)].nunique()


D = {k: load(v) for k, v in FILT.items()}

L = []
A = L.append
A('=' * 110)
A('NGC 6791 sCMD-ONLY SIX-WAY:  OG | New (dsepX_fixed) | Emul(old) | W (v1.2) | W14 (W v1.4) | W14feh (W v1.4, [Fe/H] for our Z)')
A('=' * 110)
A('')
A('Identical pipeline, identical KS cut (0.0468, 3 sigma of the original bootstrap), identical DM/E(B-V) grid.')
A('Only the isochrones differ. All legs restricted to MC sets with KTTAU 0 or 5 (W v1.4 is not trained on KTTAU = 1).')
A('')
A(f'{"leg":<11}{"rows kept":>11}{"MC sets":>10}{"best KS":>10}{"sets kept / fitted":>22}')
for k in LEGS:
    d = D[k]
    nf = n_fitted(k)
    frac = f'{d.simulation_id.nunique()}/{nf} = {d.simulation_id.nunique() / nf:.1%}' if nf == nf else '-'
    A(f'{k:<11}{len(d):>11,}{d.simulation_id.nunique():>10,}{d.ks_statistic.min():>10.4f}{frac:>22}')
A('')
A('-' * 110)
A('PARAMETER POSTERIORS')
A('-' * 110)
A(f'{"Parameter":<16}' + ''.join(f'{k:>19}' for k in LEGS))
PARS = [('Age (Gyr)', 'age_gyr', 3), ('[Fe/H]', 'FeH', 3), ('Dist. Mod.', 'distance_modulus', 3), ('E(B-V)', 'extinction_ebv', 3),
        ('He (Y)', 'He Abundance', 4), ('[alpha/Fe]', 'alphafe', 3), ('Mixing Length', 'CMIXLA', 3)]
for lab, col, f in PARS:
    row = f'{lab:<16}'
    for k in LEGS:
        d = D[k]
        x = d[col].values if col in d else IV.loc[d.simulation_id, col].values
        m = np.average(x, weights=d.w)
        sd = np.sqrt(np.average((x - m) ** 2, weights=d.w))
        row += f'{f"{m:.{f}f} +/- {sd:.{f}f}":>19}'
    A(row)
A('')
A('-' * 110)
A('AGE: paired difference on MC sets common to each pair')
A('-' * 110)
pairs = [('New', 'OG'), ('Emul(old)', 'OG'), ('W', 'OG'), ('W14', 'OG'), ('W14', 'New'), ('W14', 'W'),
         ('W14feh', 'OG'), ('W14feh', 'New'), ('W14feh', 'W14')]
A(f'{"comparison":<20}{"delta age":>12}{"SE":>9}{"z":>8}{"95% CI":>24}')
verdict = []
for a, b in pairs:
    da, db = D[a], D[b]
    common = np.intersect1d(da.simulation_id.unique(), db.simulation_id.unique())
    ca = da[da.simulation_id.isin(common)]
    cb = db[db.simulation_id.isin(common)]
    ga = ca.groupby('simulation_id').indices
    gb = cb.groupby('simulation_id').indices
    va, wa = ca.age_gyr.values, ca.w.values
    vb, wb = cb.age_gyr.values, cb.w.values
    obs = wmean(va, wa, age_bins, age_ctr) - wmean(vb, wb, age_bins, age_ctr)
    sets = np.asarray(sorted(common))
    diffs = np.empty(N_BOOT)
    for i in range(N_BOOT):
        pick = RNG.integers(0, len(sets), len(sets))
        ia = np.concatenate([ga[sets[j]] for j in pick])
        ib = np.concatenate([gb[sets[j]] for j in pick])
        diffs[i] = wmean(va[ia], wa[ia], age_bins, age_ctr) - wmean(vb[ib], wb[ib], age_bins, age_ctr)
    se = diffs.std(ddof=1)
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    z = obs / se if se > 0 else np.nan
    A(f'{f"{a} - {b}":<20}{obs:>+12.4f}{se:>9.4f}{z:>8.2f}{f"[{lo:+.4f}, {hi:+.4f}]":>24}   n={len(common)}')
    verdict.append((f'{a} - {b}', z))
A('')
A('-' * 110)
A('VERDICT (age)')
A('-' * 110)
for n, z in verdict:
    A(f'  {n:<20} z = {z:+5.2f}   {"SIGNIFICANT" if abs(z) >= 2 else "consistent with zero"}')
A('')
A('=' * 110)
t = '\n'.join(L)
open(OUT, 'w').write(t + '\n')
print(t)
print(f'\nSaved -> {OUT}')
