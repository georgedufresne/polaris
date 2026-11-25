#!/bin/bash
#SBATCH --job-name=tucmr
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=24:00:00
# One task of the replica grid: units task, task + NTASKS, ... of runs.txt x their ages (replica_worker_tuc.py).
# Submit twice (run.sh): --account=physics --array=0-239 with OFFSET=0, --account=free --array=0-31 with OFFSET=240.
# EXTRA (optional): further worker arguments, e.g. --runs/--seed/--out (see run.sh).
# Resubmitting resumes: units already in results/task_NNN.csv are skipped.
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONNOUSERSITE=1
E=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_martin_replica_mc
cd $E
/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python $E/replica_worker_tuc.py \
    --task $((SLURM_ARRAY_TASK_ID + ${OFFSET:?})) --ntasks ${NTASKS:?} ${EXTRA:-}
