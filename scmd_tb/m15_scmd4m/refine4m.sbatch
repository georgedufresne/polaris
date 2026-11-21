#!/bin/bash
# 4M-star KS descent workers (refine4m.py): each task works through the units of the first OUT until none are
# left, then the second. Arguments: OUT directories.
#SBATCH --job-name=m15_4m
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=20:00:00
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONNOUSERSITE=1 HOME=/dartfs-hpc/scratch/f007058
export TUC_AGE_LIMITS=/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
D=/dartfs-hpc/scratch/f007058/scmd_tb/m15_scmd4m
echo "task ${SLURM_ARRAY_TASK_ID} on $(hostname) start $(date '+%F %T')"
for OUT in "$@"; do
    "$PY" $D/refine4m.py "$OUT" worker --max-hours 18 2>&1 | grep -v Warn
done
echo "end $(date '+%F %T')"
