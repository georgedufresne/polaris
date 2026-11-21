#!/bin/bash
# one (variant, run) fit per array task; tasks.txt lines: "VARIANT RUN"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
export OMP_NUM_THREADS=1 HOME=/dartfs-hpc/scratch/f007058
E=/dartfs-hpc/scratch/f007058/scmd_tb/m15_errchain
read VAR RUN < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $E/tasks.txt)
cd $E/code && /dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python $E/chain.py fit $VAR $RUN 2>&1 | grep -v Warn
