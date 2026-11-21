#!/bin/bash
# Repeat of the replica grid with a different seed and 1000 runs (runs_1000.txt: the first 500 again plus 10500-10999):
# gp_minimize random_state crc32("run:age:1"); MakeTestCMD reseeds from the clock as always. Output results_s1/.
# 272 tasks (240 physics + 32 free). Rerun to resume (finished units are skipped).
set -euo pipefail
E=/dartfs-hpc/scratch/f007058/scmd_tb/m15_martin_replica
cd $E
mkdir -p logs results_s1
X="--runs $E/runs_1000.txt --seed 1 --out $E/results_s1"
J1=$(sbatch --parsable --account=physics --array=0-239 --job-name=m15mr_s1 --export=ALL,OFFSET=0,NTASKS=272,EXTRA="$X" \
     --output=$E/logs/m15mr_s1_%A_%a.out array.sbatch)
J2=$(sbatch --parsable --account=free --array=0-31 --job-name=m15mr_s1 --export=ALL,OFFSET=240,NTASKS=272,EXTRA="$X" \
     --output=$E/logs/m15mr_s1_%A_%a.out array.sbatch)
echo "$(date '+%F %T') submitted seed-1 repeat, 1000 runs (results_s1/): physics $J1, free $J2" | tee -a run.log
