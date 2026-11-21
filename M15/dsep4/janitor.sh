#!/bin/bash
# Cancel m15_d4mc array tasks that have finished (log shows "run R end") but are stuck
# in RUNNING because the spank-auks epilog crashes (as tuc_dsep4_batch/janitor.sh).
# Exits when jobs 10754215/10754216 have nothing left pending or running.
D=/dartfs-hpc/scratch/f007058/M15/dsep4
J=10754215,10754216
while squeue -j $J -h -r 2>/dev/null | grep -q . ; do
  L=""
  for id in $(squeue -j $J -h -r -t R -o "%i" 2>/dev/null); do
    j=${id%_*}; r=${id#*_}
    grep -aq "^run $r end" $D/logs/m15_d4mc_${j}_$r.out 2>/dev/null && L="$L $id"
  done
  [ -n "$L" ] && { scancel $L; echo "$(date '+%F %T') cancelled $(echo $L | wc -w): $L"; }
  sleep 600
done
echo "$(date '+%F %T') arrays finished, janitor exiting"
