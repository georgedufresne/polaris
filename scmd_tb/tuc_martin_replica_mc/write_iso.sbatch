#!/bin/bash
#SBATCH --account=free
#SBATCH --job-name=tucmr_iso
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=02:00:00
#SBATCH --output=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_martin_replica_mc/logs/write_martin_isochrones_tuc.out
cd /dartfs-hpc/scratch/f007058/scmd_tb/tuc_martin_replica_mc
export HOME=/dartfs-hpc/scratch/f007058 OMP_NUM_THREADS=1
/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python write_martin_isochrones_tuc.py
