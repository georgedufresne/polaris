#!/usr/bin/env python3
"""sCMD-ONLY age / parameter posteriors.

Derived from analyze_scmd_weighted_posteriors.py by removing every DEB leg:
  kept    - section 1 (load sCMD), 1a-bis (KS figure), 1b (initial_vars),
            1c (sCMD-only parameter posteriors), 3 (sCMD-only summary table)
  removed - sections 2, 4, 4b, 4c, 4d (DEB chi2, sCMD x DEB, DEB-only),
            5-8 (observable comparison, bar graphs, bootstrap chi2 CDF,
            overlaid age distribution) - all require DEB posteriors.

Paths are taken from env vars so the same file serves OG and new runs:
  SCMD_HOME       - dir holding filtered_ks_scmd.csv / initial_vars.csv
  SCMD_OUTPUT_DIR - where results are written
"""
import os
#!/usr/bin/env python3
"""
sCMD x DEB Weighted Age Posteriors for Observable Combinations

For all 24 variants (9 tests x 4 observable combos, with allowed_combos):
- Load sCMD data with KS 3-sigma filtering
- Load DEB chi2 per observable combo
- Merge on (mc_set, age)
- Combine weights: P_KS(mc_set) x L_DEB(age, combo)
- Generate posteriors and summary statistics

Additional parameter posteriors:
- [Fe/H]: sCMD-only, sCMD x DEB, DEB-only
- Distance modulus, E(B-V): sCMD-only, sCMD x DEB
- He abundance (Y): sCMD-only, sCMD x DEB
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.ndimage import gaussian_filter1d
from concurrent.futures import ThreadPoolExecutor
import os

# ============================================================================
# CONFIGURATION
# ============================================================================
OUTPUT_DIR = os.environ.get('SCMD_OUTPUT_DIR', '/dartfs-hpc/rc/home/8/f007058/DEB_PLOTS')
os.makedirs(OUTPUT_DIR, exist_ok=True)
HOME = os.environ.get('SCMD_HOME', '/dartfs-hpc/rc/home/8/f007058')

# Toggle for expensive plots
PLOT_KS_FULL_DISTRIBUTION = False   # two-panel KS figure (loads ~5M row CSV)

# 9 tests (5 original + 4 with V20 Secondary)
DEB_TESTS = [
    {'label': 'V20_only'},
    {'label': 'V18_pri_V20'},
    {'label': 'V18_only'},
    {'label': 'V18_sec_V20'},
    {'label': 'V18_V20'},
    {'label': 'V20_both',          'allowed_combos': ['MR']},
    {'label': 'V18_pri_V20_both',  'allowed_combos': ['MR']},
    {'label': 'V18_sec_V20_both',  'allowed_combos': ['MR']},
    {'label': 'V18_V20_all',       'allowed_combos': ['MR']},
]

# Handy sublists for plotting
ORIGINAL_TESTS = [t['label'] for t in DEB_TESTS[:5]]
V20S_TESTS = [t['label'] for t in DEB_TESTS[5:]]

# Only ML observable combination
OBSERVABLE_COMBOS = [
    {'name': 'ML',  'desc': 'M+L only'},
]

# Define age grids: 0.1 Gyr and 0.2 Gyr bin widths
AGE_BIN_CONFIGS = [
    {'width': 0.1, 'bins': np.arange(3.95, 12.15, 0.1), 'label': '0.1Gyr'},
    {'width': 0.2, 'bins': np.arange(3.85, 12.25, 0.2), 'label': '0.2Gyr'},
]
for _cfg in AGE_BIN_CONFIGS:
    _cfg['centers'] = (_cfg['bins'][:-1] + _cfg['bins'][1:]) / 2

# Default (used for stats tables etc.) — 0.1 Gyr
age_bins = AGE_BIN_CONFIGS[0]['bins']
age_centers = AGE_BIN_CONFIGS[0]['centers']

# Bins for additional parameters
feh_bins = np.arange(-0.505, 1.005, 0.02)     # [Fe/H] range
feh_centers = (feh_bins[:-1] + feh_bins[1:]) / 2
dm_bins = np.arange(10.005, 16.005, 0.02)     # distance modulus
dm_centers = (dm_bins[:-1] + dm_bins[1:]) / 2
ebv_bins = np.arange(-0.005, 0.805, 0.01)     # E(B-V) extinction
ebv_centers = (ebv_bins[:-1] + ebv_bins[1:]) / 2
he_bins = np.arange(0.195, 0.405, 0.002)      # He abundance (Y)
he_centers = (he_bins[:-1] + he_bins[1:]) / 2
afe_bins = np.array([-0.3, -0.1, 0.1, 0.3, 0.5])  # [alpha/Fe] (discrete: -0.2, 0.0, 0.2, 0.4)
afe_centers = np.array([-0.2, 0.0, 0.2, 0.4])
cmix_bins = np.arange(0.995, 2.505, 0.01)     # Mixing length (CMIXLA)
cmix_centers = (cmix_bins[:-1] + cmix_bins[1:]) / 2

# Parameter info for loop-based processing
PARAM_INFO = {
    'feh': {'bins': feh_bins, 'centers': feh_centers, 'col': 'FeH', 'label': '[Fe/H]'},
    'dm':  {'bins': dm_bins,  'centers': dm_centers,  'col': 'distance_modulus', 'label': 'Distance Modulus'},
    'ebv': {'bins': ebv_bins, 'centers': ebv_centers, 'col': 'extinction_ebv', 'label': 'E(B-V)'},
    'he':  {'bins': he_bins,  'centers': he_centers,  'col': 'He Abundance', 'label': 'He Abundance (Y)'},
    'afe': {'bins': afe_bins, 'centers': afe_centers, 'col': 'alphafe', 'label': r'[$\alpha$/Fe]', 'discrete': True},
    'cmix': {'bins': cmix_bins, 'centers': cmix_centers, 'col': 'CMIXLA', 'label': 'Mixing Length'},
}
# sCMD x DEB gets 6 params; DEB-only gets FeH, He, and Mixing Length
SCMD_DEB_PARAMS = ['feh', 'dm', 'ebv', 'he', 'afe', 'cmix']
DEB_ONLY_PARAMS = ['feh', 'he', 'cmix']

# ============================================================================
# HELPER FUNCTION: COMPUTE STATISTICS (defined early for use in loop)
# ============================================================================
def compute_stats(vals, prob, discrete=False):
    """Compute mean, median, mode, 16th/84th percentiles, and std."""
    prob = prob / np.sum(prob)  # ensure normalized
    mean = np.sum(vals * prob)
    cum = np.cumsum(prob)

    def interp_pct(cum, vals, pct):
        idx = np.searchsorted(cum, pct)
        if idx == 0: return vals[0]
        if idx >= len(vals): return vals[-1]
        if discrete:
            return vals[idx]
        frac = (pct - cum[idx-1]) / (cum[idx] - cum[idx-1] + 1e-10)
        return vals[idx-1] + frac * (vals[idx] - vals[idx-1])

    median = interp_pct(cum, vals, 0.50)
    p16 = interp_pct(cum, vals, 0.16)
    p84 = interp_pct(cum, vals, 0.84)
    mode = vals[np.argmax(prob)]
    std = np.sqrt(np.sum(prob * (vals - mean)**2))

    return {
        'mean': mean, 'median': median, 'mode': mode,
        'p16': p16, 'p84': p84, 'std': std,
    }

# ============================================================================
# 1. LOAD sCMD DATA
# ============================================================================
print("=" * 80)
print("LOADING sCMD DATA")
print("=" * 80)

# Load pre-filtered sCMD data (already cut to KS <= 0.0468 by prefilter_ks.py)
filtered_scmd = pd.read_csv(f'{HOME}/filtered_ks_scmd.csv')
bootstrap_ks = pd.read_csv(f'{HOME}/Bootstrap/bootstrap_2d_ks_results.csv')

print(f"  Loaded pre-filtered sCMD: {len(filtered_scmd):,} rows, "
      f"{filtered_scmd['simulation_id'].nunique():,} MC sets")

# Compute KS KDE (needed for probability density weighting)
ks_stats_bootstrap = bootstrap_ks['ks_statistic'].values
ks_kde = gaussian_kde(ks_stats_bootstrap, bw_method='scott')

ks_mean = np.mean(ks_stats_bootstrap)
ks_std = np.std(ks_stats_bootstrap)
ks_3sigma_limit = ks_mean + 3 * ks_std

print(f"  KS bootstrap: mean={ks_mean:.4f}, std={ks_std:.4f}, 3-sigma={ks_3sigma_limit:.4f}")

# Add KS probability density and mc_set_number
filtered_scmd['KS_Probability_Density'] = ks_kde(filtered_scmd['ks_statistic'].values)
filtered_scmd['mc_set_number'] = filtered_scmd['simulation_id']

# Keep only best (lowest KS) distance+reddening combo per (MC set, age)
n_before_best = len(filtered_scmd)
best_idx = filtered_scmd.groupby(['simulation_id', 'age_gyr'])['ks_statistic'].idxmin()
filtered_scmd = filtered_scmd.loc[best_idx]
print(f"  Best DM+EBV per (MC set, age): {len(filtered_scmd)} rows "
      f"(from {n_before_best}, {len(filtered_scmd)/n_before_best*100:.1f}%)")

# Normalize KS weights
filtered_scmd['P_KS_normalized'] = (
    filtered_scmd['KS_Probability_Density'] /
    filtered_scmd['KS_Probability_Density'].sum()
)

# ============================================================================
# 1a-bis. KS STATISTIC TWO-PANEL FIGURE (all vs good-fitting)
# ============================================================================
if PLOT_KS_FULL_DISTRIBUTION:
    print("\n  Loading full (unfiltered) KS results for two-panel plot...")
    UNFILTERED_CSV = f'{HOME}/combined_2d_ks_results_20251230_000335.csv'
    #UNFILTERED_CSV = f'{HOME}/combined_2d_ks_results_20260506_234430.csv'
    #combined_best_fits_20260506_234430
    if os.path.exists(UNFILTERED_CSV):
        ks_all = pd.read_csv(UNFILTERED_CSV, usecols=['ks_statistic'],
                             dtype={'ks_statistic': np.float32})
        ks_all_vals = ks_all['ks_statistic'].dropna().values
        ks_good_vals = filtered_scmd['ks_statistic'].values

        print(f"    All realizations: {len(ks_all_vals):,}")
        print(f"    Good-fitting (KS <= {ks_3sigma_limit:.4f}): {len(ks_good_vals):,}")

        fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(16, 6))

        # Left panel: ALL sCMDs (log x-axis to show the full range)
        ax_left.hist(ks_all_vals, bins=200, density=True, color='gray',
                     edgecolor='none', alpha=0.7)
        ax_left.axvline(ks_3sigma_limit, color='red', linestyle='--', linewidth=2.5,
                        label=f'3$\\sigma$ = {ks_3sigma_limit:.4f}')
        ax_left.set_xscale('log')
        ax_left.set_ylim(0, 10)
        ax_left.set_xlabel('KS Statistic', fontsize=18)
        ax_left.set_ylabel('Probability Density', fontsize=18)
        ax_left.tick_params(axis='both', which='major', labelsize=16, width=2.0, length=10)
        ax_left.tick_params(axis='x', which='minor', width=1.5, length=6)
        for spine in ax_left.spines.values():
            spine.set_linewidth(2.0)
        ax_left.legend(fontsize=14)
        ax_left.grid(True, alpha=0.3)

        # Right panel: good-fitting only (linear x-axis, zoomed)
        ks_bootstrap_mean = ks_mean  # bootstrap mean from earlier KDE section
        ax_right.hist(ks_good_vals, bins=80, density=True, color='steelblue',
                      edgecolor='black', linewidth=0.3, alpha=0.7)
        ax_right.axvline(ks_3sigma_limit, color='red', linestyle='--', linewidth=2.5,
                         label=f'3$\\sigma$ = {ks_3sigma_limit:.4f}')
        ax_right.axvline(ks_bootstrap_mean, color='black', linestyle='-', linewidth=2.0,
                         label=f'Bootstrap Mean = {ks_bootstrap_mean:.4f}')
        ax_right.set_xlabel('KS Statistic', fontsize=18)
        ax_right.set_ylabel('Probability Density', fontsize=18)
        ax_right.tick_params(axis='both', labelsize=16, width=2.0, length=7)
        for spine in ax_right.spines.values():
            spine.set_linewidth(2.0)
        ax_right.legend(fontsize=14)
        ax_right.grid(True, alpha=0.3)

        plt.tight_layout()
        ks_plot_path = f'{OUTPUT_DIR}/ks_statistic_two_panel.png'
        plt.savefig(ks_plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"    Saved: {ks_plot_path}")

        # --- Smooth KDE version (no bars) ---
        fig2, (ax2_left, ax2_right) = plt.subplots(1, 2, figsize=(16, 6))

        # Left panel: KDE of all KS statistics (log x-axis)
        # Use histogram to get values, then plot as filled curve
        hist_all, bin_edges_all = np.histogram(ks_all_vals, bins=500, density=True)
        bin_ctrs_all = (bin_edges_all[:-1] + bin_edges_all[1:]) / 2
        ax2_left.fill_between(bin_ctrs_all, hist_all, alpha=0.4, color='gray')
        ax2_left.plot(bin_ctrs_all, hist_all, color='black', linewidth=1.2)
        ax2_left.axvline(ks_3sigma_limit, color='red', linestyle='--', linewidth=2.5,
                         label=f'3$\\sigma$ = {ks_3sigma_limit:.4f}')
        ax2_left.set_xscale('log')
        ax2_left.set_ylim(0, 10)
        ax2_left.set_xlabel('KS Statistic', fontsize=18)
        ax2_left.set_ylabel('Probability Density', fontsize=18)
        ax2_left.tick_params(axis='both', which='major', labelsize=16, width=2.0, length=10)
        ax2_left.tick_params(axis='x', which='minor', width=1.5, length=6)
        for spine in ax2_left.spines.values():
            spine.set_linewidth(2.0)
        ax2_left.legend(fontsize=14)
        ax2_left.grid(True, alpha=0.3)

        # Right panel: KDE of good-fitting KS statistics
        kde_good = gaussian_kde(ks_good_vals, bw_method='scott')
        x_good = np.linspace(ks_good_vals.min(), ks_good_vals.max(), 500)
        y_good = kde_good(x_good)
        ax2_right.fill_between(x_good, y_good, alpha=0.4, color='steelblue')
        ax2_right.plot(x_good, y_good, color='steelblue', linewidth=2.0)
        ax2_right.axvline(ks_3sigma_limit, color='red', linestyle='--', linewidth=2.5,
                          label=f'3$\\sigma$ = {ks_3sigma_limit:.4f}')
        ax2_right.axvline(ks_bootstrap_mean, color='black', linestyle='-', linewidth=2.0,
                          label=f'Bootstrap Mean = {ks_bootstrap_mean:.4f}')
        ax2_right.set_xlabel('KS Statistic', fontsize=18)
        ax2_right.set_ylabel('Probability Density', fontsize=18)
        ax2_right.tick_params(axis='both', labelsize=16, width=2.0, length=7)
        for spine in ax2_right.spines.values():
            spine.set_linewidth(2.0)
        ax2_right.legend(fontsize=14)
        ax2_right.grid(True, alpha=0.3)

        plt.tight_layout()
        ks_plot_path2 = f'{OUTPUT_DIR}/ks_statistic_two_panel_smooth.png'
        plt.savefig(ks_plot_path2, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"    Saved: {ks_plot_path2}")

        del ks_all, ks_all_vals  # free memory
    else:
        print(f"    WARNING: {UNFILTERED_CSV} not found, skipping")

# ============================================================================
# 1b. LOAD INITIAL_VARS (all physics parameters) AND MERGE
# ============================================================================
print("\nLoading initial_vars.csv for physics parameters...")
initial_vars = pd.read_csv(f'{HOME}/initial_vars.csv')
initial_vars['file_id'] = initial_vars['file_id'].astype(int)

# Merge ALL physics columns (needed for dominance analysis later)
iv_cols = [c for c in initial_vars.columns if c != 'file_id']
filtered_scmd = pd.merge(
    filtered_scmd,
    initial_vars,
    left_on='simulation_id', right_on='file_id', how='left'
)
filtered_scmd.drop(columns=['file_id'], errors='ignore', inplace=True)

n_feh_ok = filtered_scmd['FeH'].notna().sum()
print(f"  Merged [Fe/H] + He: {n_feh_ok} / {len(filtered_scmd)} rows matched")
print(f"  FeH range [{filtered_scmd['FeH'].min():.3f}, {filtered_scmd['FeH'].max():.3f}], "
      f"He range [{filtered_scmd['He Abundance'].min():.3f}, {filtered_scmd['He Abundance'].max():.3f}]")
print(f"  DM range [{filtered_scmd['distance_modulus'].min():.3f}, {filtered_scmd['distance_modulus'].max():.3f}], "
      f"EBV range [{filtered_scmd['extinction_ebv'].min():.3f}, {filtered_scmd['extinction_ebv'].max():.3f}]")
print(f"  CMIXLA range [{filtered_scmd['CMIXLA'].min():.3f}, {filtered_scmd['CMIXLA'].max():.3f}]")

# ============================================================================
# 1c. sCMD-ONLY PARAMETER POSTERIORS (before DEB loop)
# ============================================================================
print("\n" + "=" * 80)
print("COMPUTING sCMD-ONLY PARAMETER POSTERIORS")
print("=" * 80)

scmd_only_posteriors = {}  # key: param_name -> (centers, prob_array)
scmd_only_stats = {}       # key: param_name -> stats_dict

# Age posteriors at multiple bin widths: key: (bin_label, ) or (bin_label, test, combo)
scmd_age_by_bw = {}        # key: bin_label -> (centers, prob)
scmd_age_stats_by_bw = {}  # key: bin_label -> stats_dict

scmd_weights = filtered_scmd['P_KS_normalized'].values

# sCMD-only: Age (compute for all bin widths)
hist_age_scmd, _ = np.histogram(filtered_scmd['age_gyr'].values, bins=age_bins, weights=scmd_weights)
if hist_age_scmd.sum() > 0:
    hist_age_scmd = hist_age_scmd / hist_age_scmd.sum()
    scmd_only_posteriors['age'] = (age_centers, hist_age_scmd)
    scmd_only_stats['age'] = compute_stats(age_centers, hist_age_scmd)
    print(f"  Age:  mean={scmd_only_stats['age']['mean']:.2f}, mode={scmd_only_stats['age']['mode']:.2f} Gyr")

for cfg in AGE_BIN_CONFIGS:
    h, _ = np.histogram(filtered_scmd['age_gyr'].values, bins=cfg['bins'], weights=scmd_weights)
    if h.sum() > 0:
        h = h / h.sum()
        scmd_age_by_bw[cfg['label']] = (cfg['centers'], h)
        scmd_age_stats_by_bw[cfg['label']] = compute_stats(cfg['centers'], h)

# sCMD-only: [Fe/H]
feh_vals = filtered_scmd['FeH'].dropna()
feh_w = filtered_scmd.loc[feh_vals.index, 'P_KS_normalized'].values
hist_feh_scmd, _ = np.histogram(feh_vals.values, bins=feh_bins, weights=feh_w)
if hist_feh_scmd.sum() > 0:
    hist_feh_scmd = hist_feh_scmd / hist_feh_scmd.sum()
    scmd_only_posteriors['feh'] = (feh_centers, hist_feh_scmd)
    scmd_only_stats['feh'] = compute_stats(feh_centers, hist_feh_scmd)
    print(f"  [Fe/H]: mean={scmd_only_stats['feh']['mean']:.3f}, mode={scmd_only_stats['feh']['mode']:.3f}")

# sCMD-only: Distance Modulus
hist_dm_scmd, _ = np.histogram(filtered_scmd['distance_modulus'].values, bins=dm_bins, weights=scmd_weights)
if hist_dm_scmd.sum() > 0:
    hist_dm_scmd = hist_dm_scmd / hist_dm_scmd.sum()
    scmd_only_posteriors['dm'] = (dm_centers, hist_dm_scmd)
    scmd_only_stats['dm'] = compute_stats(dm_centers, hist_dm_scmd)
    print(f"  DM:   mean={scmd_only_stats['dm']['mean']:.3f}, mode={scmd_only_stats['dm']['mode']:.3f}")

# sCMD-only: E(B-V)
hist_ebv_scmd, _ = np.histogram(filtered_scmd['extinction_ebv'].values, bins=ebv_bins, weights=scmd_weights)
if hist_ebv_scmd.sum() > 0:
    hist_ebv_scmd = hist_ebv_scmd / hist_ebv_scmd.sum()
    scmd_only_posteriors['ebv'] = (ebv_centers, hist_ebv_scmd)
    scmd_only_stats['ebv'] = compute_stats(ebv_centers, hist_ebv_scmd)
    print(f"  E(B-V): mean={scmd_only_stats['ebv']['mean']:.3f}, mode={scmd_only_stats['ebv']['mode']:.3f}")

# sCMD-only: He Abundance
he_vals = filtered_scmd['He Abundance'].dropna()
he_w = filtered_scmd.loc[he_vals.index, 'P_KS_normalized'].values
hist_he_scmd, _ = np.histogram(he_vals.values, bins=he_bins, weights=he_w)
if hist_he_scmd.sum() > 0:
    hist_he_scmd = hist_he_scmd / hist_he_scmd.sum()
    scmd_only_posteriors['he'] = (he_centers, hist_he_scmd)
    scmd_only_stats['he'] = compute_stats(he_centers, hist_he_scmd)
    print(f"  He(Y): mean={scmd_only_stats['he']['mean']:.4f}, mode={scmd_only_stats['he']['mode']:.4f}")

# sCMD-only: [alpha/Fe]
afe_vals = filtered_scmd['alphafe'].dropna()
afe_w = filtered_scmd.loc[afe_vals.index, 'P_KS_normalized'].values
hist_afe_scmd, _ = np.histogram(afe_vals.values, bins=afe_bins, weights=afe_w)
if hist_afe_scmd.sum() > 0:
    hist_afe_scmd = hist_afe_scmd / hist_afe_scmd.sum()
    scmd_only_posteriors['afe'] = (afe_centers, hist_afe_scmd)
    scmd_only_stats['afe'] = compute_stats(afe_centers, hist_afe_scmd, discrete=True)
    print(f"  [a/Fe]: mean={scmd_only_stats['afe']['mean']:.3f}, mode={scmd_only_stats['afe']['mode']:.3f}")

# sCMD-only: Mixing Length
cmix_vals = filtered_scmd['CMIXLA'].dropna()
cmix_w = filtered_scmd.loc[cmix_vals.index, 'P_KS_normalized'].values
hist_cmix_scmd, _ = np.histogram(cmix_vals.values, bins=cmix_bins, weights=cmix_w)
if hist_cmix_scmd.sum() > 0:
    hist_cmix_scmd = hist_cmix_scmd / hist_cmix_scmd.sum()
    scmd_only_posteriors['cmix'] = (cmix_centers, hist_cmix_scmd)
    scmd_only_stats['cmix'] = compute_stats(cmix_centers, hist_cmix_scmd)
    print(f"  Mixing Length: mean={scmd_only_stats['cmix']['mean']:.4f}, mode={scmd_only_stats['cmix']['mode']:.4f}")

# ============================================================================
# 3. sCMD-ONLY SUMMARY TABLE
# ============================================================================
print("\n\n" + "=" * 80)
print("sCMD-ONLY PARAMETER SUMMARY")
print("=" * 80)

scmd_only_summary_data = []
for pname in ['age', 'feh', 'dm', 'ebv', 'he', 'afe', 'cmix']:
    if pname in scmd_only_stats:
        s = scmd_only_stats[pname]
        lbl = {'age': 'Age (Gyr)', 'feh': '[Fe/H]', 'dm': 'Dist. Mod.',
               'ebv': 'E(B-V)', 'he': 'He (Y)', 'afe': r'[$\alpha$/Fe]', 'cmix': 'Mixing Length'}[pname]
        scmd_only_summary_data.append({
            'Parameter': lbl, 'Mean': s['mean'], 'Median': s['median'],
            'Mode': s['mode'], 'Std': s['std'], 'P16': s['p16'], 'P84': s['p84'],
        })
        print(f"  {lbl:>12s}: mean={s['mean']:.4f}, median={s['median']:.4f}, "
              f"mode={s['mode']:.4f}, std={s['std']:.4f}")

scmd_only_summary_df = pd.DataFrame(scmd_only_summary_data)
csv_scmd_only_params = f'{OUTPUT_DIR}/scmd_only_parameter_summary.csv'
scmd_only_summary_df.to_csv(csv_scmd_only_params, index=False)
print(f"\n  Saved: {csv_scmd_only_params}")