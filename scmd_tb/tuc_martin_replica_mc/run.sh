#!/bin/bash
# Submit the replica grid: NP physics (240 minus my running physics CPUs) + 32 free tasks, each working through every
# (NP + 32)th (run, age) unit.
# Writes the job ids to run.log. Rerun to resume after a time-out (finished units are skipped).
set -euo pipefail
E=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_martin_replica_mc
cd $E
mkdir -p logs results
NP=$((240 - $(squeue -u f007058 -h -t R -A physics -o "%C" | awk "{s+=\$1} END {print s+0}")))
J1=$(sbatch --parsable --account=physics --array=0-$((NP - 1)) --export=ALL,OFFSET=0,NTASKS=$((NP + 32)) \
     --output=$E/logs/tucmr_%A_%a.out array.sbatch)
J2=$(sbatch --parsable --account=free --array=0-31 --export=ALL,OFFSET=$NP,NTASKS=$((NP + 32)) \
     --output=$E/logs/tucmr_%A_%a.out array.sbatch)
echo "$(date '+%F %T') submitted replica grid: physics $J1, free $J2" | tee -a run.log
