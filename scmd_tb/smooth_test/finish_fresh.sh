#!/bin/bash
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export HOME=/dartfs-hpc/scratch/f007058
until [ "$(ls /dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/newmodel_fresh/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 232 ]; do sleep 20; done
sleep 20; scancel 10198095 10198096 2>/dev/null
cd /dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/newmodel_fresh/code && $PY -m isochronetoolbox.scmd bootstrap-collect --out /dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/newmodel_fresh 2>&1 | grep -v Warn
