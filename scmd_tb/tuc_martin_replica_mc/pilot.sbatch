#!/bin/bash
#SBATCH --job-name=tucmr_pilot
#SBATCH --account=physics
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_martin_replica_mc/logs/pilot_%j.out
# Pilot: three units through the full replica (MakeTestCMD -> verbatim KS -> gp_minimize), for timing and memory.
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONNOUSERSITE=1
E=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_martin_replica_mc
cd $E
/usr/bin/time -v /dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python $E/replica_worker_tuc.py --task 0 --ntasks 1 \
    --units 10567:10200,10567:10400,10105:12200 --out $E/pilot_results
