#!/bin/bash
# Waits for the 47 Tuc DSEP4 Monte Carlo grid to finish (finish_tail.sh writes runs_mcfix_built.txt), then archives
# its iso files and isochrones (archive_iso.sh tuc4). PID in after_tuc4.pid.
B=/dartfs-hpc/scratch/f007058/tuc_dsep4_batch; A=/dartfs-hpc/rc/lab/C/ChaboyerB/george/polaris/archive
until [ -s $B/runs_mcfix_built.txt ]; do sleep 900; done
sleep 600   # let the last isochrone builds settle
cd $A && echo "$(date '+%F %T') 47 Tuc grid done; $(sbatch archive_iso.sh tuc4)" >> $A/after_tuc4.log
