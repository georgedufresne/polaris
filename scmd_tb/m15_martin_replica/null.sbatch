#!/bin/bash
#SBATCH --job-name=m15mr_null
#SBATCH --account=physics
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=03:00:00
# Bootstrap null of the replica (null_worker.py): 50 tasks x 200 fake data sets from the best unit.
# Submit: sbatch --array=0-49 --export=ALL,RUN=<run>,AGE=<Myr>,DM=<dm> --output=logs/null_%A_%a.out null.sbatch
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONNOUSERSITE=1
E=/dartfs-hpc/scratch/f007058/scmd_tb/m15_martin_replica
cd $E
/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python $E/null_worker.py --task $SLURM_ARRAY_TASK_ID --n 200 \
    --run ${RUN:?} --age ${AGE:?} --dm ${DM:?} ${EXTRA:-}
