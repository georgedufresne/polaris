#!/bin/bash
# Wait for the fixed DSEP4 rerun (M15/dsep4/run_m15_dsep4mc.batch, arrays in M15/dsep4/jobids_mcfix.txt) to finish,
# then: cancel tasks that hang after logging "end" (spank-auks crash at exit), list the runs with isochrones, compute
# their age limits, and run this experiment's run.sh.
D=/dartfs-hpc/scratch/f007058/M15/dsep4
O=/dartfs-hpc/scratch/f007058/scmd_tb/m15_final_mc
I=/dartfs-hpc/scratch/f007058/M15/iso_dsep4mc
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
echo "waiting for the DSEP4 rerun ($(date '+%F %T'))"
ended() { grep -al "^run [0-9]* end" $D/logs_mcfix/m15_d4fix_*_*.out 2>/dev/null | wc -l; }
while true; do
  n=$(ended); q=$(squeue -u f007058 -h -n m15_d4fix | wc -l)
  [ "$n" -ge 500 ] && break
  [ "$q" -eq 0 ] && { echo "queue empty with $n of 500 runs ended"; break; }
  sleep 600
done
echo "DSEP4 rerun: $(ended) of 500 runs ended ($(date '+%F %T')); jobs left: $(squeue -u f007058 -h -n m15_d4fix | wc -l)"
scancel -u f007058 -n m15_d4fix 2>/dev/null; sleep 30
ls -d $I/1*/ 2>/dev/null | while read d; do r=$(basename $d); [ -s $I/$r/isos/isochrones_example.txt ] && echo $r; done | sort -n > $O/runs.txt
echo "runs with isochrones: $(wc -l < $O/runs.txt) (list in runs.txt)"
$PY $O/age_limits_m15mc.py $O/runs.txt $O/m15_age_limits_d4mc.json || exit 1
exec bash $O/run.sh
