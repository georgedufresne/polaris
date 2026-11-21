#!/bin/bash
# Every 10 min: cancel W14 array tasks that are still RUNNING but whose log already has the pipeline's
# "done" line (Slurm spank-auks hang at task exit). Stops when the array leaves the queue. PID in reap_zombies.pid.
J=11106027; P=/dartfs-hpc/scratch/f007058/scmd_W14_pipeline
while squeue -h -j $J 2>/dev/null | grep -q .; do
  ids=""
  for t in $(squeue -j $J -h -t R -o "%K"); do
    grep -aqE "^done |already done" $P/logs/scmd_${J}_$t.out 2>/dev/null && ids="$ids ${J}_$t"
  done
  [ -n "$ids" ] && scancel $ids && echo "$(date '+%F %T') cancelled finished zombies:$ids" >> $P/zombies_cancelled.txt
  sleep 600
done
echo "$(date '+%F %T') array gone, reaper exits" >> $P/zombies_cancelled.txt
