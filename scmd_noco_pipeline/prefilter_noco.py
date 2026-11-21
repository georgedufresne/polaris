#!/usr/bin/env python3
"""
Pre-filter the NOCO combined_2d_ks_results by the KS 3-sigma cut.

Repointed copy of ~/prefilter_ks.py.  The bootstrap null is NOT recomputed:
the user directed us to reuse the existing one, so KS_3SIGMA is the same
0.0468 (bootstrap mean=0.0268, std=0.0067) used for the original data.
Keeping the threshold fixed is also what makes OG vs NEW comparable.

Input  : newest combined_2d_ks_results_*.csv in /dartfs-hpc/scratch/f007058/analysis_noco
Output : <that dir>/filtered_ks_scmd.csv
"""

import glob
import os

import numpy as np
import pandas as pd

ANALYSIS_DIR = os.environ.get('NOCO_DIR', '/dartfs-hpc/scratch/f007058/analysis_noco')

cands = sorted(glob.glob(os.path.join(ANALYSIS_DIR, 'combined_2d_ks_results_*.csv')),
               key=os.path.getmtime)
if not cands:
    raise SystemExit(f"No combined_2d_ks_results_*.csv in {ANALYSIS_DIR} -- run merge_noco.py first")
INPUT_CSV = cands[-1]
OUTPUT_CSV = os.path.join(ANALYSIS_DIR, 'filtered_ks_scmd.csv')

KS_3SIGMA = 0.0468

USE_COLS = ['ks_statistic', 'simulation_id', 'age_gyr',
            'distance_modulus', 'extinction_ebv']
DTYPES = {
    'ks_statistic': np.float32,
    'simulation_id': np.int32,
    'age_gyr': np.float32,
    'distance_modulus': np.float32,
    'extinction_ebv': np.float32,
}

print(f"Reading {INPUT_CSV} (columns: {USE_COLS})...")
df = pd.read_csv(INPUT_CSV, usecols=USE_COLS, dtype=DTYPES)
print(f"  Loaded {len(df):,} rows, {df['simulation_id'].nunique():,} MC sets")

mask = df['ks_statistic'] <= KS_3SIGMA
df_filtered = df[mask].copy()
print(f"  After KS <= {KS_3SIGMA} cut: {len(df_filtered):,} rows "
      f"({len(df_filtered)/len(df)*100:.1f}%)")
print(f"  Unique MC sets surviving: {df_filtered['simulation_id'].nunique():,}")

df_filtered.to_csv(OUTPUT_CSV, index=False)
print(f"  Saved {OUTPUT_CSV} ({os.path.getsize(OUTPUT_CSV)/1e6:.1f} MB)")
