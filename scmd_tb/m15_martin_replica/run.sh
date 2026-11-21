#!/bin/bash
# Submit the replica grid: 272 tasks (240 physics + 32 free), each working through every 272nd (run, age) unit.
# Writes the job ids to run.log. Rerun to resume after a time-out (finished units are skipped).
set -euo pipefail
E=/dartfs-hpc/scratch/f007058/scmd_tb/m15_martin_replica
cd $E
mkdir -p logs results
J1=$(sbatch --parsable --account=physics --array=0-239 --export=ALL,OFFSET=0,NTASKS=272 \
     --output=$E/logs/m15mr_%A_%a.out array.sbatch)
J2=$(sbatch --parsable --account=free --array=0-31 --export=ALL,OFFSET=240,NTASKS=272 \
     --output=$E/logs/m15mr_%A_%a.out array.sbatch)
echo "$(date '+%F %T') submitted replica grid: physics $J1, free $J2" | tee -a run.log
