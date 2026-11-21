#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=01:30:00
export OMP_NUM_THREADS=1 HOME=/dartfs-hpc/scratch/f007058
G=/dartfs-hpc/scratch/f007058/scmd_tb/m15_grid_d4
RUN=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $G/escan/runs.txt)
cd $G/code && /dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python $G/escan/escan_fine.py $G $RUN
