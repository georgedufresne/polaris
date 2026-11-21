#!/bin/bash
# All 1000 47 Tuc DSEP4 runs with the Monte Carlo read fixed (run_tuc_dsep4_mcfix.batch), within 240 physics + 32 free:
#  - while the M15 production refit runs (driver pid in scmd_tb/m15_final_mc/driver.pid, 200 physics + 32 free), this array
#    gets 20 physics and the NGC 6791 W14 sCMD array (11106027) 20;
#  - after it: this array 200 physics, W14 40, and the 32 runs held on the free account (10968-10999) are released;
#  - when the 20 runs of d3_d4_20_panels_hrd.png have isochrones, the figure is redrawn (scmd_tb/tuc_iso_compare/mcfix);
#  - at the end: hung tasks cancelled (they log "end" first), runs list written, W14 given the whole physics share (240).
B=/dartfs-hpc/scratch/f007058/tuc_dsep4_batch
I=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso_d4v2_mcfix
C=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso_compare
M15D=$(cat /dartfs-hpc/scratch/f007058/scmd_tb/m15_final_mc/driver.pid)
W14=11106027
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
echo $$ > $B/run_mcfix1000.pid
mkdir -p $B/logs_mcfix/mass /dartfs-hpc/scratch/f007058/Tracks_tuc_dsep4_mcfix $I
log() { echo "$(date '+%F %T') $*" | tee -a $B/jobids_mcfix.txt; }
J1=$(sbatch --parsable --account=physics --array=10000-10967%20 $B/run_tuc_dsep4_mcfix.batch)
J2=$(sbatch --parsable --hold --account=free --array=10968-10999 $B/run_tuc_dsep4_mcfix.batch)
log "47 Tuc 1000-run fix: physics $J1 (10000-10967, 20 at once while M15 refits), free $J2 (10968-10999, held until M15 is done); M15 driver $M15D"
nend() { grep -al "^run [0-9]* end" $B/logs_mcfix/tuc_d4fix_*_*.out 2>/dev/null | wc -l; }
FIG="10000 10002 $(seq 10004 10021)"
nfig() { n=0; for r in $FIG; do [ -s $I/$r/isos/isochrones_example.txt ] && n=$((n+1)); done; echo $n; }
phase=1; plotted=0
while true; do
  if [ $phase = 1 ] && ! kill -0 $M15D 2>/dev/null; then
    sleep 120
    scontrol update JobId=$J1 ArrayTaskThrottle=200 >/dev/null 2>&1
    scontrol update JobId=$W14 ArrayTaskThrottle=40 >/dev/null 2>&1
    scontrol release $J2
    phase=2; log "M15 refit finished: $J1 -> 200 physics, W14 -> 40, free $J2 released"
  fi
  if [ $plotted = 0 ] && [ "$(nfig)" -ge 20 ]; then
    (cd $C && $PY plot_20_panels_mcfix.py) && plotted=1 && log "redrawn 20-panel figures in $C/mcfix"
  fi
  [ "$(nend)" -ge 1000 ] && break
  squeue -h -j $J1,$J2 2>/dev/null | grep -q . || break
  sleep 600
done
scancel $J1 $J2 2>/dev/null
ls -d $I/1*/ 2>/dev/null | while read d; do r=$(basename $d); [ -s $I/$r/isos/isochrones_example.txt ] && echo $r; done | sort -n > $B/runs_mcfix_built.txt
log "47 Tuc fix done: $(nend) runs ended, $(wc -l < $B/runs_mcfix_built.txt) with isochrones (runs_mcfix_built.txt)"
[ $plotted = 0 ] && (cd $C && $PY plot_20_panels_mcfix.py)
scontrol update JobId=$W14 ArrayTaskThrottle=240 >/dev/null 2>&1 && log "W14 array $W14 -> 240"
