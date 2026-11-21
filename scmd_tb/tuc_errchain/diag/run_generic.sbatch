#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:40:00
export OMP_NUM_THREADS=2 HOME=/dartfs-hpc/scratch/f007058
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
B=/dartfs-hpc/scratch/f007058/scmd_tb; T=$B/tuc_errchain
cd $T/diag


$PY diag_generic.py tuc.json 2>&1 | grep -v Warn
