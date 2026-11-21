#!/bin/bash
# Wait for the ../m15_final_yingiso driver to exit (it holds all 272 CPUs), then run this experiment's run.sh.
P=$(cat /dartfs-hpc/scratch/f007058/scmd_tb/m15_final_yingiso/driver.pid)
echo "waiting for m15_final_yingiso driver $P ($(date '+%F %T'))"
while kill -0 $P 2>/dev/null; do sleep 120; done
sleep 120
echo "m15_final_yingiso driver exited ($(date '+%F %T')); jobs left in queue: $(squeue -u f007058 -h | wc -l)"
exec bash /dartfs-hpc/scratch/f007058/scmd_tb/m15_final_phxbc/run.sh
