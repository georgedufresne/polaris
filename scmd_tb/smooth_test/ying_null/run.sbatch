#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=02:00:00
export OMP_NUM_THREADS=1 PYTHONNOUSERSITE=1
$SCMD_PYTHON /dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/ying_null/ying_null.py /dartfs-hpc/scratch/f007058/scmd_tb/smooth_dev ${SLURM_ARRAY_TASK_ID} 31
