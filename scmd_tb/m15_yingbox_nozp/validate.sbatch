#!/bin/bash
#SBATCH --job-name=yb_validate
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
export OMP_NUM_THREADS=1 PYTHONNOUSERSITE=1 HOME=/dartfs-hpc/scratch/f007058
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
B=/dartfs-hpc/scratch/f007058/scmd_tb
TUC_AGE_LIMITS=/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json $PY $B/m15_yingbox_nozp/validate_yingbox.py m15 2>&1 | grep -v Warn
TUC_AGE_LIMITS=$B/tuc_final/tuc_age_limits_d4.json $PY $B/tuc_yingbox_nozp/validate_yingbox.py tuc 2>&1 | grep -v Warn
