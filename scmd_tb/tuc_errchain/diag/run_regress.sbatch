#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=00:40:00
export OMP_NUM_THREADS=2 HOME=/dartfs-hpc/scratch/f007058
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
B=/dartfs-hpc/scratch/f007058/scmd_tb; T=$B/tuc_errchain
cd $T/diag
$PY regress_ridge.py $B/tuc_grid_d4/code $B/tuc_grid_d4/setup.json $T/diag/reg_prod.npy 2>&1 | grep -v Warn
$PY regress_ridge.py $T/code $T/setup.json $T/diag/reg_exp.npy ridge 2>&1 | grep -v Warn
$PY -c "import numpy as np; a=np.load('$T/diag/reg_prod.npy'); b=np.load('$T/diag/reg_exp.npy'); print('REGRESSION identical:', a.shape==b.shape and np.array_equal(a,b))"
