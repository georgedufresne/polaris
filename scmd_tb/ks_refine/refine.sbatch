#!/bin/bash
# KS refinement workers (refine.py): each task refines units of the first grid until none are left, then
# of the second. Arguments: pairs of OUT AGE_LIMITS_JSON, optionally followed by FIRST_RUNS files via env
# REFINE_FIRST_1 / REFINE_FIRST_2.
#SBATCH --job-name=ks_refine
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=20:00:00
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONNOUSERSITE=1 HOME=/dartfs-hpc/scratch/f007058
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
D=/dartfs-hpc/scratch/f007058/scmd_tb/ks_refine
echo "task ${SLURM_ARRAY_TASK_ID} on $(hostname) start $(date '+%F %T')"
n=1
while [ $# -ge 2 ]; do
    OUT=$1; LIM=$2; shift 2
    F=REFINE_FIRST_$n; FIRST=${!F}
    "$PY" $D/refine.py "$OUT" worker --age-limits "$LIM" --max-hours 12 ${FIRST:+--first $FIRST} 2>&1 | grep -v Warn
    n=$((n + 1))
done
echo "end $(date '+%F %T')"
