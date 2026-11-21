#!/bin/bash
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export HOME=/dartfs-hpc/scratch/f007058
until [ "$(ls /dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/newmodel/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 200 ] && [ -f /dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/newmodel/data_vs_model.csv ]; do sleep 30; done
sleep 20; scancel 10197773 10197774 2>/dev/null
cd /dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/newmodel/code && $PY -m isochronetoolbox.scmd bootstrap-collect --out /dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/newmodel 2>&1 | grep -v Warn
$PY -c "
import pandas as pd
d = pd.read_csv('/dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/newmodel/data_vs_model.csv'); print(f'data vs new model, {len(d)} realisations: KS {d.ks_statistic.mean():.4f} +- {d.ks_statistic.std():.4f}')"
