#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=12G
#SBATCH --time=01:00:00
export OMP_NUM_THREADS=1 HOME=/dartfs-hpc/scratch/f007058
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
B=/dartfs-hpc/scratch/f007058/scmd_tb; D=$B/literature/impact
echo "=== M15"
cd $B/m15_grid_d4/code && $PY $D/impact.py $B/m15_grid_d4 \
  '[[10286,14.1,15.48,0.08],[10387,13.3,15.49,0.08],[10145,14.4,15.48,0.08],[10392,13.3,15.52,0.065]]' \
  '[[-0.45,0.02],[-1.0,0.02],[-0.45,0.017],[-1.0,0.017]]' $D/m15_impact.csv 2>&1 | grep -v Warn
