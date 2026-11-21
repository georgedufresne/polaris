#!/bin/bash
# Tail of the 47 Tuc 1000-run fix. Some runs hit the 24 h limit (with Monte Carlo on, masses near 0.8-0.95 Msun can take
# 6-13 h). run_tuc_dsep4_mcfix.batch skips masses whose track already exists, so a resubmission finishes only what is
# missing. Every 10 min:
#  - a run with no "run R end" line that is neither running nor pending is resubmitted (physics, --time=72:00:00),
#    at most twice;
#  - tasks that logged "end" but hang in RUNNING are cancelled;
#  - the NGC 6791 W14 array (11106027) gets ArrayTaskThrottle = 240 - (running 47 Tuc tasks).
# When all 1000 runs have ended: runs_mcfix_built.txt is written and W14 gets 240.
B=/dartfs-hpc/scratch/f007058/tuc_dsep4_batch
I=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso_d4v2_mcfix
W14=11106027
echo $$ > $B/finish_tail.pid
declare -A tries
log() { echo "$(date '+%F %T') $*" | tee -a $B/jobids_mcfix.txt; }
ended() { grep -ahoE "^run [0-9]+ end" $B/logs_mcfix/tuc_d4fix_*_*.out 2>/dev/null | awk '{print $2}' | sort -un; }
last=-1
while true; do
  E=$(ended); ne=$(echo "$E" | grep -c .)
  # hung tasks (logged end, still RUNNING)
  for t in $(squeue -u f007058 -h -n tuc_d4fix -t R -o "%i"); do r=${t#*_}; echo "$E" | grep -qx "$r" && scancel $t; done
  [ "$ne" -ge 1000 ] && break
  Q=$(squeue -u f007058 -h -r -n tuc_d4fix -o "%i" | sed 's/.*_//' | sort -un)
  todo=""
  for r in $(seq 10000 10999); do
    echo "$E" | grep -qx "$r" && continue
    echo "$Q" | grep -qx "$r" && continue
    n=${tries[$r]:-0}; [ "$n" -ge 2 ] && continue
    tries[$r]=$((n + 1)); todo="$todo,$r"
  done
  if [ -n "$todo" ]; then
    J=$(sbatch --parsable --account=physics --time=72:00:00 --array=${todo#,} $B/run_tuc_dsep4_mcfix.batch)
    log "resubmitted (72 h) runs ${todo#,} as $J"
  fi
  nr=$(squeue -u f007058 -h -n tuc_d4fix -t R,PD -r | wc -l)
  thr=$((240 - nr)); [ "$thr" -lt 20 ] && thr=20
  if [ "$thr" != "$last" ]; then scontrol update JobId=$W14 ArrayTaskThrottle=$thr >/dev/null 2>&1; last=$thr; log "W14 throttle $thr (47 Tuc tasks queued/running: $nr; runs ended: $ne)"; fi
  sleep 600
done
ls -d $I/1*/ 2>/dev/null | while read d; do r=$(basename $d); [ -s $I/$r/isos/isochrones_example.txt ] && echo $r; done | sort -n > $B/runs_mcfix_built.txt
scontrol update JobId=$W14 ArrayTaskThrottle=240 >/dev/null 2>&1
log "47 Tuc fix complete: 1000 runs ended, $(wc -l < $B/runs_mcfix_built.txt) with isochrones; W14 throttle 240"
