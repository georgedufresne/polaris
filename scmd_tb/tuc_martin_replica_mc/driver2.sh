#!/bin/bash
# Second driver (2026-11-26), after the first grid (try1_sparse/) failed on our sparse EEP isochrones.
# 0. write the densely resampled isochrones (write_iso.sbatch) and wait;
# 1. grid: run.sh; after 400 finished units, stop everything if more than 2% of sCMDs have < 1,000,000 stars (the
#    try1 symptom); otherwise wait for all units of ages_per_run.json, resubmitting once (run.sh resumes) if tasks are
#    gone with units missing. Tasks that logged "finished" but stay RUNNING (spank-auks) are cancelled;
# 2. analyze_replica_tuc.py without null -> logs/analyze_grid.out; 3. null at the best unit; 4. analysis and
#    compare_pipelines.txt. Stops (exit 1) instead of carrying on with missing units or a missing null.
# Log: driver2.log. PID: driver.pid.
set -u
E=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_martin_replica_mc
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
cd $E
echo $$ > driver.pid
reap() {
  squeue -h -u f007058 -n "$1" -t R -o "%F %K" | while read -r A T; do
    f=$(ls logs/*_${A}_${T}.out 2>/dev/null | head -1)
    [ -n "$f" ] && grep -a -q "finished in\|^task .*: base" "$f" && scancel "${A}_${T}" && echo "  reaped ${A}_${T} ($(date '+%F %T'))"
  done
}
done_units() { cat results/task_*.csv 2>/dev/null | grep -v -c '^simulation_id' ; }
small_frac() { $PY -c "
import glob, pandas as pd
f = glob.glob('results/task_*.csv')
g = pd.concat([pd.read_csv(x) for x in f]) if f else pd.DataFrame({'n_scmd': []})
print(len(g), round(float((g.n_scmd < 1e6).mean()) if len(g) else 0.0, 4))"; }

echo "=== isochrones $(date '+%F %T')"
J=$(sbatch --parsable write_iso.sbatch); echo "write_iso job $J"
while squeue -h -j $J 2>/dev/null | grep -q .; do sleep 30; done
cat logs/write_martin_isochrones_tuc.out
grep -q "age blocks in" logs/write_martin_isochrones_tuc.out || { echo "isochrone writer failed"; exit 1; }
NEXP=$($PY -c "import json; d=json.load(open('ages_per_run.json')); r=[int(x) for x in open('runs.txt').read().split()]; print(sum(len(d[str(k)]) for k in r))")

echo "=== grid $(date '+%F %T'): $NEXP units"
./run.sh
checked=0; resub=0
while :; do
  sleep 300
  n=$(done_units)
  if [ $checked -eq 0 ] && [ "$n" -ge 400 ]; then
    read N F < <(small_frac)
    echo "  check after $N units: fraction of sCMDs below 1e6 stars = $F ($(date '+%F %T'))"
    if $PY -c "import sys; sys.exit(0 if float('$F') > 0.02 else 1)"; then
      scancel -n tucmr -u f007058; echo "STOP: too many small sCMDs"; exit 1
    fi
    checked=1
  fi
  [ "$n" -ge "$NEXP" ] && break
  if ! squeue -h -u f007058 -n tucmr | grep -q .; then
    if [ $resub -eq 0 ]; then echo "  tasks gone with $n of $NEXP units: resubmitting ($(date '+%F %T'))"; ./run.sh; resub=1
    else echo "STOP: tasks gone again with $n of $NEXP units"; exit 1; fi
  fi
  reap tucmr
done
sleep 60; reap tucmr
echo "grid done $(date '+%F %T'): $(done_units) of $NEXP units; $(small_frac)"

echo "=== analysis without null $(date '+%F %T')"
$PY analyze_replica_tuc.py > logs/analyze_grid.out 2>&1 || true
head -3 logs/analyze_grid.out
read RUN AGE DM < <($PY -c "import pandas as pd; g=pd.read_csv('combined_results.csv'); b=g.loc[g.ks_statistic.idxmin()]; print(int(b.simulation_id), int(b.age_myr), f'{b.distance_modulus:.6f}')")
echo "=== null at run $RUN, $AGE Myr, DM $DM $(date '+%F %T')"
mkdir -p null
NJ=$(sbatch --parsable --array=0-49 --export=ALL,RUN=$RUN,AGE=$AGE,DM=$DM --output=$E/logs/null_%A_%a.out null.sbatch)
echo "null job $NJ" | tee -a run.log
while :; do
  sleep 120
  [ "$(ls null/null_task_*.csv 2>/dev/null | wc -l)" -ge 50 ] && break
  squeue -h -u f007058 -n tucmr_null | grep -q . || break
  reap tucmr_null
done
sleep 30; reap tucmr_null
NN=$(ls null/null_task_*.csv 2>/dev/null | wc -l)
echo "null done $(date '+%F %T'): $NN of 50 parts"
[ "$NN" -ge 45 ] || { echo "STOP: null incomplete; see logs/null_${NJ}_*.out"; exit 1; }

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
