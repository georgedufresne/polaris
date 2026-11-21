#!/bin/bash
#SBATCH --job-name=d4v2_iso
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=02:00:00
# task i builds lines i+1, i+101, i+201, ... of runs_complete.txt (build_tuc_d4_iso_v2.sh)
D=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_d4v2_build
export HOME=/dartfs-hpc/scratch/f007058 OMP_NUM_THREADS=1
for r in $(awk -v i=$SLURM_ARRAY_TASK_ID '(NR - 1) % 100 == i' $D/runs_complete.txt); do
  $D/build_tuc_d4_iso_v2.sh $r
done
echo "task $SLURM_ARRAY_TASK_ID done $(date '+%F %T')"
