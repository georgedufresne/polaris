#!/bin/bash
# Bootstrap null of the 2D KS statistic, one slice of the realisations per array
# task (bootstrap.py).  Point it at a finished grid directory (its setup.json
# fixes the model); the cell defaults to the grid's lowest KS.  Give every array
# the total number of slices and its offset, so arrays on two accounts never
# take the same slice:
#
#   sbatch --account=physics --array=0-239 --output=$OUT/logs/boot_%A_%a.out \
#          scripts/scmd_bootstrap.sbatch $OUT --n 10000 --ntasks 272 --task-offset 0
#   sbatch --account=free    --array=0-31  --output=$OUT/logs/boot_%A_%a.out \
#          scripts/scmd_bootstrap.sbatch $OUT --n 10000 --ntasks 272 --task-offset 240
#   python -m isochronetoolbox.scmd bootstrap-collect --out $OUT
#
# Kerberos: as scmd_worker.sbatch -- with SCMD_PYTHON set and everything on scratch.
#
#SBATCH --job-name=scmd_boot
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=12:00:00

OUT=$1
shift

if [ -n "${SCMD_PYTHON}" ]; then
    PY=${SCMD_PYTHON}
else
    module load python/anaconda3 2>/dev/null
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate "${SCMD_ENV:-fidanka-env}"
    PY=python
fi
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONNOUSERSITE=1

cd "${SLURM_SUBMIT_DIR}"
echo "task ${SLURM_ARRAY_TASK_ID} on $(hostname)  start $(date '+%F %T')"
"${PY}" -m isochronetoolbox.scmd bootstrap --out "$OUT" "$@"
echo "end $(date '+%F %T')"
