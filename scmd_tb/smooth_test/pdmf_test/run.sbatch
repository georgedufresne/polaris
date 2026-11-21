#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=01:00:00
export OMP_NUM_THREADS=1 PYTHONNOUSERSITE=1
$SCMD_PYTHON /dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/pdmf_test/mini.py /dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/pdmf_test ${SLURM_ARRAY_TASK_ID}
