#!/bin/bash
# After 278 KTTAU=1 tasks of the W14 sCMD array were cancelled (kt1_cancelled.txt), run_W14.sh can no longer reach its
# 1000 "done" logs. When every remaining task (1000 - 278 = 722) has logged done, cancel the array (tasks can hang in
# RUNNING after exiting) so run_W14.sh leaves its wait loop and runs the analysis. Does not touch run_W14.sh.
P=/dartfs-hpc/scratch/f007058/scmd_W14_pipeline
J=11106027
NEED=$((1000 - $(tr ',' '\n' < $P/kt1_cancelled.txt | grep -c .)))
echo "waiting for $NEED done logs of $J ($(date '+%F %T'))"
ndone() { grep -alE "^done |already done" $P/logs/scmd_${J}_*.out 2>/dev/null | wc -l; }
while [ "$(ndone)" -lt "$NEED" ]; do
  squeue -h -j $J 2>/dev/null | grep -q . || { echo "array left the queue at $(ndone) done"; exit 0; }
  sleep 900
done
scancel $J 2>/dev/null
echo "all $NEED remaining tasks done; array cancelled so run_W14.sh continues ($(date '+%F %T'))"
