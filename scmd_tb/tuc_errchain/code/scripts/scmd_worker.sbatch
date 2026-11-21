#!/bin/bash
# One sCMD worker per array task.  Workers claim units of work from a shared
# run list (atomic mkdir), so several arrays -- on different accounts -- can
# be pointed at the same OUT and run list without duplicating anything.
#
# Submit from the repository root (or a copy of it); size each array to the
# CPUs you may use:
#
#   OUT=/dartfs-hpc/scratch/$USER/scmd_out; mkdir -p $OUT/logs
#   sbatch --account=physics --array=0-239 --output=$OUT/logs/%A_%a.out \
#          scripts/scmd_worker.sbatch $OUT runs.txt
#   sbatch --account=free    --array=0-31  --output=$OUT/logs/%A_%a.out \
#          scripts/scmd_worker.sbatch $OUT runs.txt
#
# Extra arguments after the run list go to the worker
# (e.g. --exact-age, --bc data/ngc6791/bc_mist_gaia_edr3_fehp025.npz).
#
# Kerberos: on Discovery /dartfs-hpc/rc (home, lab shares) is mounted
# sec=krb5 and batch jobs get their ticket from auks, which fails when many
# tasks start at once (and expires after ~10 h); /dartfs-hpc/scratch is
# sec=sys.  So this script is plain bash (a login shell would read ~/.bashrc)
# and, with SCMD_PYTHON set, runs that interpreter directly -- keep it, the
# code, the isochrones (--iso) and OUT on scratch and no task needs a ticket.
#
#SBATCH --job-name=scmd
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00

OUT=$1
RUNS=$2
shift 2

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
"${PY}" -m isochronetoolbox.scmd worker --out "$OUT" --runs "$RUNS" --max-hours 8 "$@"
echo "end $(date '+%F %T')"
