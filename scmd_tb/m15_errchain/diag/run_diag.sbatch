#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=01:00:00
export OMP_NUM_THREADS=2 HOME=/dartfs-hpc/scratch/f007058
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
B=/dartfs-hpc/scratch/f007058/scmd_tb; D=$B/m15_errchain/diag
cd $D
$PY regress.py $B/m15_grid_d4/code $B/m15_grid_d4/setup.json $D/reg_prod.npy 2>&1 | grep -v Warn
$PY regress.py $B/m15_errchain/code $B/m15_errchain/setup.json $D/reg_exp.npy 2>&1 | grep -v Warn
$PY -c "import numpy as np; a=np.load('$D/reg_prod.npy'); b=np.load('$D/reg_exp.npy'); print('REGRESSION identical:', a.shape==b.shape and np.array_equal(a,b))"
$PY diag.py 2>&1 | grep -v Warn
