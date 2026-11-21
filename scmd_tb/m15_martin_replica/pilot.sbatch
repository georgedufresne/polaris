#!/bin/bash
#SBATCH --job-name=m15mr_pilot
#SBATCH --account=physics
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:30:00
#SBATCH --output=/dartfs-hpc/scratch/f007058/scmd_tb/m15_martin_replica/logs/pilot_%j.out
# Pilot: three units through the full replica (MakeTestCMD -> verbatim KS -> gp_minimize), for timing and memory.
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONNOUSERSITE=1
E=/dartfs-hpc/scratch/f007058/scmd_tb/m15_martin_replica
cd $E
/usr/bin/time -v /dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python $E/replica_worker.py --task 0 --ntasks 1 \
    --units 10089:12600,10089:13400,10027:12800 --out $E/pilot_results
