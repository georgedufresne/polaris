#!/bin/bash
#SBATCH --job-name=v4_fine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=03:00:00
# one (run, age) per task: line SLURM_ARRAY_TASK_ID + 1 of cells.txt
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONNOUSERSITE=1 HOME=/dartfs-hpc/scratch/f007058
D=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_v4_fine
read RUN AGE <<< "$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $D/cells.txt)"
cd $D
/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python fine_grid.py $RUN $AGE 2>&1 | grep -v -i warn
