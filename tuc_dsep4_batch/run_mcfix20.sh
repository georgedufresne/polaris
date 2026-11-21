#!/bin/bash
# The 20 runs of d3_d4_20_panels_hrd.png with DSEP4 reading its Monte Carlo inputs, then the redrawn figure.
# CPU limit 240 physics: lower the NGC 6791 W14 array to 20 concurrent, hold these 20 runs and release them once my
# running physics CPUs are <= 220; restore W14 to 40 at the end.
B=/dartfs-hpc/scratch/f007058/tuc_dsep4_batch
C=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso_compare
I=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso_d4v2_mcfix
W14=11106027
echo $$ > $B/run_mcfix20.pid
mkdir -p $B/logs_mcfix/mass
scontrol update JobId=$W14 ArrayTaskThrottle=20 && echo "W14 array $W14 throttled to 20 $(date '+%F %T')"
J=$(sbatch --parsable --hold --account=physics --array=10000,10002,10004-10021%20 $B/run_tuc_dsep4_mcfix.batch)
echo "tuc_d4fix array $J held $(date '+%F %T')" | tee -a $B/jobids_mcfix.txt
phys() { squeue -u f007058 -h -A physics -t R -o '%C' | awk '{s+=$1} END{print s+0}'; }
until [ "$(phys)" -le 220 ]; do sleep 300; done
scontrol release $J && echo "released at $(phys) physics CPUs $(date '+%F %T')"
nend() { grep -al "^run [0-9]* end" $B/logs_mcfix/tuc_d4fix_${J}_*.out 2>/dev/null | wc -l; }
while [ "$(nend)" -lt 20 ]; do
  squeue -h -j $J 2>/dev/null | grep -q . || break
  sleep 600
done
scancel $J 2>/dev/null
echo "runs ended: $(nend)/20; isochrones: $(ls $I/*/isos/isochrones_example.txt 2>/dev/null | wc -l) $(date '+%F %T')"
scontrol update JobId=$W14 ArrayTaskThrottle=40 && echo "W14 throttle back to 40"
cd $C && /dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python plot_20_panels_mcfix.py
echo "all done $(date '+%F %T')"
