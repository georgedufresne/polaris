#!/usr/bin/env python3
"""
Write a side-by-side text comparison of the sCMD-only posteriors from the
ORIGINAL grid and the NEW (noco / fixed-physics) grid.

Both summaries are produced by analyze_scmd_only.py; the same KS 3-sigma
threshold (0.0468, from the original bootstrap) was applied to both, so the
numbers are directly comparable.
"""
import os
import pandas as pd

OG   = os.environ.get('OG_DIR',   '/dartfs-hpc/scratch/f007058/analysis_OG')
NEW  = os.environ.get('NEW_DIR',  '/dartfs-hpc/scratch/f007058/analysis_noco')
OUT  = os.environ.get('CMP_OUT',  '/dartfs-hpc/scratch/f007058/scmd_og_vs_noco_comparison.txt')

og  = pd.read_csv(f'{OG}/scmd_only_parameter_summary.csv').set_index('Parameter')
new = pd.read_csv(f'{NEW}/scmd_only_parameter_summary.csv').set_index('Parameter')


def counts(d):
    """rows surviving the KS cut and the number of distinct MC sets."""
    f = pd.read_csv(f'{d}/filtered_ks_scmd.csv', usecols=['simulation_id'])
    return len(f), f['simulation_id'].nunique()


try:
    og_rows, og_sets = counts(OG)
except Exception:
    og_rows, og_sets = counts('/dartfs-hpc/rc/home/8/f007058')
new_rows, new_sets = counts(NEW)

L = []
L.append('=' * 92)
L.append('sCMD-ONLY PARAMETER POSTERIORS  --  ORIGINAL GRID vs NEW (fixed-physics, no overshoot)')
L.append('=' * 92)
L.append('')
L.append(f'  original : {OG}')
L.append(f'  new      : {NEW}')
L.append(f'  produced by analyze_scmd_only.py (no DEB leg, no bootstrap recomputed)')
L.append(f'  KS cut   : ks_statistic <= 0.0468  (3 sigma of the original bootstrap null,')
L.append(f'             mean 0.0268 / std 0.0067) -- identical for both, so they compare directly')
L.append('')
L.append(f'  rows surviving KS cut :  OG {og_rows:>9,}      NEW {new_rows:>9,}')
L.append(f'  distinct MC sets      :  OG {og_sets:>9,}      NEW {new_sets:>9,}')
L.append('')

for stat in ['Mean', 'Median', 'Mode', 'Std']:
    L.append('-' * 92)
    L.append(f'{stat.upper()}')
    L.append('-' * 92)
    L.append(f'{"Parameter":<18}{"original":>14}{"new":>14}{"new - orig":>14}{"shift / orig sigma":>22}')
    for p in og.index:
        if p not in new.index:
            continue
        a, b = og.loc[p, stat], new.loc[p, stat]
        d = b - a
        sig = og.loc[p, 'Std']
        rel = f'{d / sig:+.2f} sigma' if sig > 0 else 'n/a'
        L.append(f'{p:<18}{a:>14.4f}{b:>14.4f}{d:>+14.4f}{rel:>22}')
    L.append('')

L.append('-' * 92)
L.append('68% CREDIBLE INTERVAL  (P16 - P84)')
L.append('-' * 92)
L.append(f'{"Parameter":<18}{"original":>26}{"new":>26}{"width chg":>14}')
for p in og.index:
    if p not in new.index:
        continue
    oi = f'[{og.loc[p, "P16"]:.4f}, {og.loc[p, "P84"]:.4f}]'
    ni = f'[{new.loc[p, "P16"]:.4f}, {new.loc[p, "P84"]:.4f}]'
    dw = ((new.loc[p, 'P84'] - new.loc[p, 'P16']) -
          (og.loc[p, 'P84'] - og.loc[p, 'P16']))
    L.append(f'{p:<18}{oi:>26}{ni:>26}{dw:>+14.4f}')
L.append('')
L.append('=' * 92)

txt = '\n'.join(L)
with open(OUT, 'w') as fh:
    fh.write(txt + '\n')
print(txt)
print(f'\nSaved -> {OUT}')
