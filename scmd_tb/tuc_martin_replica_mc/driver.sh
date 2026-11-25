#!/bin/bash
# Driver for the 47 Tuc replica of M. Ying's 2D-KS pipeline on our Monte Carlo-on isochrones (2026-11-25).
# 1. grid: run.sh (NP physics + 32 free tasks of array.sbatch -> replica_worker_tuc.py), wait until every unit of
#    ages_per_run.json is in results/ or no tucmr task is left. Tasks that logged "finished" but stay RUNNING (Slurm
#    spank-auks crash at exit, seen in ../m15_martin_replica) are cancelled one by one.
# 2. analyze_replica_tuc.py without a null (best KS by age) -> logs/analyze_grid.out.
# 3. null: null.sbatch, 50 tasks x 200 sets at the best unit (base fit DM), -> null/.
# 4. analyze_replica_tuc.py with the null -> logs/analyze_replica_tuc.out; compare_pipelines.txt.
# Log: driver.log. PID: driver.pid.
set -u
E=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_martin_replica_mc
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
cd $E
echo $$ > driver.pid
NEXP=$($PY -c "import json; d=json.load(open('ages_per_run.json')); r=[int(x) for x in open('runs.txt').read().split()]; print(sum(len(d[str(k)]) for k in r))")
done_units() { cat results/task_*.csv 2>/dev/null | grep -v -c '^simulation_id' ; }
reap() {   # cancel array tasks of job name $1 whose log already says "finished"
  squeue -h -u f007058 -n "$1" -t R -o "%F %K" | while read -r A T; do
    f=$(ls logs/*_${A}_${T}.out 2>/dev/null | head -1)
    [ -n "$f" ] && grep -a -q "finished\|^task .*: base" "$f" && scancel "${A}_${T}" && echo "  reaped ${A}_${T} ($(date '+%F %T'))"
  done
}

echo "=== grid $(date '+%F %T'): $NEXP units"
./run.sh
sleep 120
while :; do
  n=$(done_units)
  squeue -h -u f007058 -n tucmr | grep -q . || break
  [ "$n" -ge "$NEXP" ] && break
  reap tucmr
  sleep 300
done
sleep 60; reap tucmr
echo "grid done $(date '+%F %T'): $(done_units) of $NEXP units"
squeue -h -u f007058 -n tucmr | grep -q . && echo "  still queued/running tucmr tasks:" && squeue -h -u f007058 -n tucmr -o "%i %T"
[ "$(done_units)" -lt "$NEXP" ] && echo "  WARNING: units missing; rerun ./run.sh to resume (finished units are skipped)"

echo "=== analysis without null $(date '+%F %T')"
$PY analyze_replica_tuc.py > logs/analyze_grid.out 2>&1 || true
head -5 logs/analyze_grid.out
read RUN AGE DM < <($PY -c "import pandas as pd; g=pd.read_csv('combined_results.csv'); b=g.loc[g.ks_statistic.idxmin()]; print(int(b.simulation_id), int(b.age_myr), f'{b.distance_modulus:.6f}')")
echo "=== null at run $RUN, $AGE Myr, DM $DM $(date '+%F %T')"
mkdir -p null
NJ=$(sbatch --parsable --array=0-49 --export=ALL,RUN=$RUN,AGE=$AGE,DM=$DM --output=$E/logs/null_%A_%a.out null.sbatch)
echo "null job $NJ" | tee -a run.log
sleep 120
while squeue -h -u f007058 -n tucmr_null | grep -q .; do
  [ "$(ls null/null_task_*.csv 2>/dev/null | wc -l)" -ge 50 ] && break
  reap tucmr_null
  sleep 120
done
sleep 30; reap tucmr_null
echo "null done $(date '+%F %T'): $(ls null/null_task_*.csv 2>/dev/null | wc -l) of 50 parts"

echo "=== analysis $(date '+%F %T')"
$PY analyze_replica_tuc.py > logs/analyze_replica_tuc.out 2>&1
cat logs/analyze_replica_tuc.out
{
  echo "47 Tuc, Monte Carlo-on DSEP4 isochrones: M. Ying's pipeline (this folder) vs ours (../tuc_final_mc_refined)"
  echo "--- this replica (logs/analyze_replica_tuc.out)"
  grep -a -E "units \(|^null:|kde |normal |survival |fit probability" logs/analyze_replica_tuc.out | grep -v calibration
  echo "--- our pipeline (../tuc_final_mc/compare_production.txt)"
  grep tuc_final_mc_refined ../tuc_final_mc/compare_production.txt
  echo "--- Ying et al. 2025 Table 7, 47 Tuc 2D KS: 11.76 +/- 0.76 Gyr, DM 13.40, E 0.02"
} > compare_pipelines.txt
cat compare_pipelines.txt
echo "all done $(date '+%F %T')"
