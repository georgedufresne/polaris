#!/bin/bash
# Driver (2026-11-22): NGC 6791 W v1.4 leg with [Fe/H] converted to our Z convention (W14feh), then re-age every leg.
#  1. isochrones: ../emul_W14/gen_W14feh.batch (40 tasks x 4 CPUs, gen_W14_iso.py --z-factor 1.121) for the 668 KTTAU 0/5
#     runs (../emul_W14/runs_W14feh.txt); one retry for runs left without an isochrone
#  2. sCMDs: run_scmd_W14feh.batch over the runs with isochrones; array throttle = 240 physics CPUs minus my other
#     physics jobs, rechecked every 10 min; tasks that logged "done" but hang (spank-auks) are cancelled
#  3. waits for the DSEP3 (11131955) and W14 (11131934) reruns
#  4. completeness, every leg's KTTAU 0/5 runs: a run needs as many ages in its 2D-KS table as its isochrone has
#     (81); incomplete tables are moved (not deleted) to <ksdir>_partial/ and those runs rerun once; runs still
#     incomplete stay out and are listed in incomplete_<leg>.txt
#  5. merge + prefilter for New (analysis_noco), W14 (analysis_W14), W14feh (analysis_W14feh); each leg's previous
#     filtered_ks_scmd.csv is kept as filtered_ks_scmd_pre_rerun.csv
#  6. scmd_6way_compare.py (-> ../scmd_6way_comparison.txt), age_breakdown_W14feh.py, ../scmd_W14_pipeline/age_breakdown.py
# Log: run_W14feh.log. PID: driver.pid. ../scmd_tb/tuc_final_mc/run.sh waits for "all done" here.
set -u
S=/dartfs-hpc/scratch/f007058; E=$S/emul_W14; P=$S/scmd_W14feh_pipeline; PY=$S/envs/fidanka-env/bin/python
NOCO_J=11131955; W14R_J=11131934; CAP=240
cd $P; echo $$ > $P/driver.pid
log() { echo "$(date '+%F %T') $*"; }
RL=$E/runs_W14feh.txt; N=$(wc -l < $RL)

built_list() { while read r; do f=$S/Isochrones_W14feh/$r/isos/isochrones_example.txt
               [ -s $f ] && [ "$(stat -c %s $f)" -gt 1000000 ] && echo $r; done < $RL; }
# running physics CPUs of mine, excluding array $1
other_physics() { squeue -u f007058 -h -t R -A physics -o '%i %C' | grep -v "^${1}_" | awk '{s+=$2} END {print s+0}'; }
# cancel running tasks of array $1 whose log ($2 with IDX for the task index) already says done
reap() { local j=$1 pat=$2 t i f
  for t in $(squeue -h -j $j -t R -o '%i' 2>/dev/null); do i=${t#*_}; f=${pat//IDX/$i}
    grep -aqE '^done |already done' $f 2>/dev/null && { scancel $t; log "cancelled finished zombie $t" >> $P/zombies_cancelled.txt; }
  done; }
wait_array() { local j=$1 pat=$2
  while squeue -h -j $j 2>/dev/null | grep -q .; do
    reap $j "$pat"; reap $NOCO_J "$S/scmd_noco_pipeline/logs/rerun_${NOCO_J}_IDX.out"
    reap $W14R_J "$S/scmd_W14_pipeline/logs/scmd_rerun_${W14R_J}_IDX.out"
    if [ "${3:-}" = throttle ]; then
      t=$((CAP - $(other_physics $j))); [ $t -lt 1 ] && t=1
      scontrol update JobId=$j ArrayTaskThrottle=$t > /dev/null 2>&1
    fi
    sleep 600
  done; }
# incomplete runs of a leg: $1 = 2D-KS dir, $2 = run list, $3 = isochrone path pattern with RUN
incomplete() { $PY - "$1" "$2" "$3" <<'PYEOF'
import os, subprocess, sys
import pandas as pd
d, rl, iso = sys.argv[1:4]
for r in open(rl).read().split():
    f = f"{d}/2d_ks_results_{r}.csv"
    n = 0
    if os.path.exists(f) and os.path.getsize(f) > 0:
        try:
            n = pd.read_csv(f, usecols=["age_gyr"]).age_gyr.round(1).nunique()
        except Exception:
            n = 0
    if n >= 81:
        continue
    out = subprocess.run(f"awk '!/^#/ {{print $2}}' {iso.replace('RUN', r)} | sort -u | wc -l", shell=True,
                         capture_output=True, text=True).stdout.strip()
    expect = min(int(out or 81), 81)
    if n < expect:
        print(r)
PYEOF
}
move_partial() { local d=$1 r; mkdir -p ${d}_partial
  for r in $(cat $2); do for f in $d/2d_ks_results_$r.csv $d/best_fit_$r.csv; do [ -e $f ] && mv -f $f ${d}_partial/; done; done; }

log "=== 1. isochrones ($N runs)"
B=$(built_list | wc -l)
if [ "$B" -lt "$N" ]; then
  K=$(( (CAP - $(other_physics none)) / 4 )); [ $K -gt 40 ] && K=40; [ $K -lt 5 ] && K=5
  G=$(sbatch --parsable --array=0-39%$K --export=ALL,NT=40 $E/gen_W14feh.batch); log "gen array $G (%$K)" | tee -a $P/jobids.txt
  while squeue -h -j $G 2>/dev/null | grep -q .; do sleep 120; done
  B=$(built_list | wc -l); log "isochrones built: $B/$N"
  if [ "$B" -lt "$N" ]; then
    comm -23 <(sort $RL) <(built_list | sort) | sort -n > $E/runs_W14feh_retry.txt
    n2=$(wc -l < $E/runs_W14feh_retry.txt); NT2=$(( n2 < 20 ? n2 : 20 ))
    G2=$(sbatch --parsable --array=0-$((NT2 - 1)) --export=ALL,NT=$NT2,RUNLIST=runs_W14feh_retry.txt $E/gen_W14feh.batch)
    log "gen retry $G2 ($n2 runs)" | tee -a $P/jobids.txt
    while squeue -h -j $G2 2>/dev/null | grep -q .; do sleep 120; done
  fi
fi
built_list > $P/runs.txt; NB=$(wc -l < $P/runs.txt); log "isochrones built: $NB/$N (runs.txt)"
[ "$NB" -gt 0 ] || { log "no isochrones; stopping"; exit 1; }

log "=== 2. sCMDs"
T=$((CAP - $(other_physics none))); [ $T -lt 1 ] && T=1; [ $T -gt $NB ] && T=$NB
J=$(sbatch --parsable --array=0-$((NB - 1))%$T --export=ALL,RUNLIST=runs.txt $P/run_scmd_W14feh.batch)
log "sCMD array $J ($NB runs, throttle $T, adjusted to $CAP minus my other physics jobs)" | tee -a $P/jobids.txt
wait_array $J "$P/logs/scmd_${J}_IDX.out" throttle
log "sCMD array done: $(ls $S/2d_ks_W14feh/2d_ks_results_*.csv 2>/dev/null | wc -l) KS tables"

log "=== 3. waiting for the DSEP3 and W14 reruns"
while squeue -h -j $NOCO_J,$W14R_J 2>/dev/null | grep -q .; do
  reap $NOCO_J "$S/scmd_noco_pipeline/logs/rerun_${NOCO_J}_IDX.out"; reap $W14R_J "$S/scmd_W14_pipeline/logs/scmd_rerun_${W14R_J}_IDX.out"
  sleep 600
done

log "=== 4. completeness"
$PY - <<PYEOF
import pandas as pd
iv = pd.read_csv("$S/analysis_W/initial_vars.csv").set_index("file_id")
for src, out in (("$S/scmd_noco_pipeline/runs.txt", "$P/runs_noco_k05.txt"), ("$S/scmd_W14_pipeline/runs_W14.txt", "$P/runs_W14_k05.txt")):
    rr = [int(x) for x in open(src).read().split()]
    open(out, "w").write("\n".join(str(r) for r in rr if iv.KTTAU.get(r) in (0, 5)) + "\n")
PYEOF
incomplete $S/2d_ks_W14feh $P/runs.txt "$S/Isochrones_W14feh/RUN/isos/isochrones_example.txt" > $P/retry_W14feh.txt
incomplete $S/2d_ks_noco $P/runs_noco_k05.txt "$S/Isochrones_noco/RUN/isos/isochrones_example.txt" > $P/retry_noco.txt
incomplete $S/2d_ks_W14 $P/runs_W14_k05.txt "$S/Isochrones_W14/RUN/isos/isochrones_example.txt" > $P/retry_W14.txt
log "incomplete: W14feh $(wc -l < $P/retry_W14feh.txt), New $(wc -l < $P/retry_noco.txt), W14 $(wc -l < $P/retry_W14.txt)"
RJ=""
if [ -s $P/retry_W14feh.txt ]; then
  move_partial $S/2d_ks_W14feh $P/retry_W14feh.txt; cp $P/retry_W14feh.txt $P/runs_retry.txt
  j=$(sbatch --parsable --array=0-$(( $(wc -l < $P/runs_retry.txt) - 1 )) --export=ALL,RUNLIST=runs_retry.txt $P/run_scmd_W14feh.batch)
  log "W14feh retry $j" | tee -a $P/jobids.txt; RJ="$RJ $j:$P/logs/scmd_${j}_IDX.out"
fi
if [ -s $P/retry_noco.txt ]; then
  move_partial $S/2d_ks_noco $P/retry_noco.txt; cp $P/retry_noco.txt $S/scmd_noco_pipeline/runs_noco_retry.txt
  j=$(sbatch --parsable --array=0-$(( $(wc -l < $P/retry_noco.txt) - 1 )) --export=ALL,RUNLIST=runs_noco_retry.txt $S/scmd_noco_pipeline/run_scmd_noco_rerun.batch)
  log "New retry $j" | tee -a $P/jobids.txt; RJ="$RJ $j:$S/scmd_noco_pipeline/logs/rerun_${j}_IDX.out"
fi
if [ -s $P/retry_W14.txt ]; then
  move_partial $S/2d_ks_W14 $P/retry_W14.txt; cp $P/retry_W14.txt $S/scmd_W14_pipeline/runs_W14_retry.txt
  j=$(sbatch --parsable --array=0-$(( $(wc -l < $P/retry_W14.txt) - 1 )) --export=ALL,RUNLIST=runs_W14_retry.txt $S/scmd_W14_pipeline/run_scmd_W14_rerun.batch)
  log "W14 retry $j" | tee -a $P/jobids.txt; RJ="$RJ $j:$S/scmd_W14_pipeline/logs/scmd_rerun_${j}_IDX.out"
fi
for x in $RJ; do wait_array ${x%%:*} "${x#*:}"; done
incomplete $S/2d_ks_W14feh $P/runs.txt "$S/Isochrones_W14feh/RUN/isos/isochrones_example.txt" > $P/incomplete_W14feh.txt
incomplete $S/2d_ks_noco $P/runs_noco_k05.txt "$S/Isochrones_noco/RUN/isos/isochrones_example.txt" > $P/incomplete_noco.txt
incomplete $S/2d_ks_W14 $P/runs_W14_k05.txt "$S/Isochrones_W14/RUN/isos/isochrones_example.txt" > $P/incomplete_W14.txt
for leg in W14feh noco W14; do
  case $leg in W14feh) d=$S/2d_ks_W14feh;; noco) d=$S/2d_ks_noco;; W14) d=$S/2d_ks_W14;; esac
  [ -s $P/incomplete_$leg.txt ] && move_partial $d $P/incomplete_$leg.txt
done
log "still incomplete after the retry (kept out of the merge): W14feh $(wc -l < $P/incomplete_W14feh.txt), New $(wc -l < $P/incomplete_noco.txt), W14 $(wc -l < $P/incomplete_W14.txt)"

log "=== 5. merge + prefilter"
for a in analysis_noco analysis_W14; do
  [ -e $S/$a/filtered_ks_scmd_pre_rerun.csv ] || cp -p $S/$a/filtered_ks_scmd.csv $S/$a/filtered_ks_scmd_pre_rerun.csv
done
mkdir -p $S/analysis_W14feh; cp -p $S/analysis_W/initial_vars.csv $S/analysis_W14feh/
( cd $S/scmd_noco_pipeline && $PY merge_noco.py > $P/logs/merge_noco.out 2>&1 && $PY prefilter_noco.py ) || log "New merge/prefilter FAILED"
( cd $S/scmd_W14_pipeline && $PY merge_W14.py > $P/logs/merge_W14.out 2>&1 && $PY prefilter_W14.py ) || log "W14 merge/prefilter FAILED"
( cd $P && $PY merge_W14feh.py > $P/logs/merge_W14feh.out 2>&1 && $PY prefilter_W14feh.py ) || log "W14feh merge/prefilter FAILED"

log "=== 6. comparison"
[ -e $S/scmd_W14_pipeline/age_breakdown_pre_rerun.txt ] || cp -p $S/scmd_W14_pipeline/age_breakdown.txt $S/scmd_W14_pipeline/age_breakdown_pre_rerun.txt
$PY $P/scmd_6way_compare.py > $P/logs/6way.out 2>&1 || log "6-way comparison FAILED (logs/6way.out)"
$PY $P/age_breakdown_W14feh.py > $P/logs/age_breakdown_W14feh.out 2>&1 || log "W14feh age breakdown FAILED"
$PY $S/scmd_W14_pipeline/age_breakdown.py > $P/logs/age_breakdown_W14.out 2>&1 || log "W14 age breakdown FAILED (reads home; rerun by hand)"
tail -30 $S/scmd_6way_comparison.txt 2>/dev/null
log "all done"
