#!/bin/bash
# bootstrap null of one variant at its best cell; nullvars.txt lines: VARIANT
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=03:00:00
export OMP_NUM_THREADS=1 HOME=/dartfs-hpc/scratch/f007058
E=/dartfs-hpc/scratch/f007058/scmd_tb/m15_errchain
VAR=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $E/nullvars.txt)
cd $E/code && /dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python $E/chain.py null $VAR 240 2>&1 | grep -v Warn
