#!/bin/bash
# Cancel tuc_d4mc array tasks that have finished (log shows "run R end") but are
# stuck in RUNNING because the spank-auks epilog crashes (malloc_consolidate).
# Exits when jobs 10656804/10656805/10671878 have nothing left pending or running.
B=/dartfs-hpc/scratch/f007058/tuc_dsep4_batch
while squeue -j 10656804,10656805,10671878 -h -r 2>/dev/null | grep -q . ; do
  L=""
  for id in $(squeue -j 10656804,10656805,10671878 -h -r -t R -o "%i" 2>/dev/null); do
    j=${id%_*}; r=${id#*_}
    grep -aq "^run $r end" $B/logs/tuc_d4mc_${j}_$r.out 2>/dev/null && L="$L $id"
  done
  [ -n "$L" ] && { scancel $L; echo "$(date '+%F %T') cancelled $(echo $L | wc -w): $L"; }
  sleep 600
done
echo "$(date '+%F %T') arrays finished, janitor exiting"
