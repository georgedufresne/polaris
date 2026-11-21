#!/bin/bash
# Wait for each bootstrap's 272 slices, cancel tasks hung in the auks epilog, collect; then the 47 Tuc posterior.
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export HOME=/dartfs-hpc/scratch/f007058
N=/dartfs-hpc/scratch/f007058/scmd_tb/boot_ngc6791
T=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02
until [ "$(ls $N/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 272 ]; do sleep 30; done
sleep 30; scancel 10090122 10090123 2>/dev/null
cd $N/code && $PY -m isochronetoolbox.scmd bootstrap-collect --out $N 2>&1 | grep -v Warn
until [ "$(ls $T/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 272 ]; do sleep 30; done
sleep 30; scancel 10091784 10091785 2>/dev/null
cd $T/code_boot && $PY -m isochronetoolbox.scmd bootstrap-collect --out $T 2>&1 | grep -v Warn
$PY -m isochronetoolbox.scmd posterior --out $T 2>&1 | grep -v Warn
