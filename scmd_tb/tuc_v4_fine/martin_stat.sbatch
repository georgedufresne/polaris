#!/bin/bash
#SBATCH --job-name=v4_mstat
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=01:00:00
# task 0-49: null realisations 20*i .. 20*i+19;  task 50: the fits
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONNOUSERSITE=1 HOME=/dartfs-hpc/scratch/f007058
cd /dartfs-hpc/scratch/f007058/scmd_tb/tuc_v4_fine
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
if [ "$SLURM_ARRAY_TASK_ID" -eq 50 ]; then
    $PY martin_stat.py fit 2>&1 | grep -v -i warn
else
    $PY martin_stat.py null $((SLURM_ARRAY_TASK_ID * 20)) 20 2>&1 | grep -v -i warn
fi
