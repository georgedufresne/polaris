#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=01:00:00
export OMP_NUM_THREADS=1 PYTHONNOUSERSITE=1
$SCMD_PYTHON /dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/data_vs_model.py /dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/smooth_resample 10099:11.4:13.38:0.015 12
