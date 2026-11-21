#!/bin/bash
#SBATCH --job-name=yi_validate
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
export OMP_NUM_THREADS=1 PYTHONNOUSERSITE=1 HOME=/dartfs-hpc/scratch/f007058
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
$PY /dartfs-hpc/scratch/f007058/scmd_tb/m15_yingiso_nozp/validate_yingiso.py 2>&1 | grep -v Warn
