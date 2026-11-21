#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
export OMP_NUM_THREADS=2 HOME=/dartfs-hpc/scratch/f007058
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
B=/dartfs-hpc/scratch/f007058/scmd_tb; D=$B/m15_errchain/diag
cd $D
$PY 
$PY 
$PY diag2.py 2>&1 | grep -v Warn
