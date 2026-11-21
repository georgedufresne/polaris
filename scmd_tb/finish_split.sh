#!/bin/bash
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export HOME=/dartfs-hpc/scratch/f007058
until [ "$(ls /dartfs-hpc/scratch/f007058/scmd_tb/hugs/run10005_split/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 40 ] && [ "$(ls /dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02_split/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 232 ]; do sleep 30; done
sleep 30; scancel 10095043 10095044 10095045 2>/dev/null
cd /dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02_split/code && $PY -m isochronetoolbox.scmd bootstrap-collect --out /dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02_split 2>&1 | grep -v Warn
$PY /dartfs-hpc/scratch/f007058/scmd_tb/tuc_statistics.py /dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02_split 2>&1 | grep -v Warn
cd /dartfs-hpc/scratch/f007058/scmd_tb/hugs/run10005_split/code && $PY -m isochronetoolbox.scmd bootstrap-collect --out /dartfs-hpc/scratch/f007058/scmd_tb/hugs/run10005_split 2>&1 | grep -v Warn
$PY - <<'PYEOF'
import pandas as pd
H = '/dartfs-hpc/scratch/f007058/scmd_tb/hugs/run10005_split'
b = pd.read_csv(f'{H}/bootstrap_2d_ks_results.csv'); g = pd.read_csv(f'{H}/combined_2d_ks_results.csv')
v = b.p95_diff; lim = v.mean() + 3 * v.std(ddof=0)
print(f'HUGS 10005, 95th percentile: null {v.mean():.4f} +- {v.std(ddof=0):.4f}, 3 sigma {lim:.4f}; best fit {g.p95_diff.min():.4f}; rows passing {(g.p95_diff <= lim).sum()} of {len(g)}')
PYEOF
