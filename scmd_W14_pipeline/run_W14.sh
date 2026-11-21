#!/bin/bash
# NGC 6791, W v1.4 leg (2026-11-18). After the W v1.4 isochrones (../emul_W14, gen_W14 array in ../emul_W14/jobids.txt):
#  1. list the runs with complete isochrones (runs_W14.txt);
#  2. the sCMD pipeline (run_scmd_W14.batch = ../scmd_emul_pipeline/run_scmd_W_physics.batch with the W14 folders) on the
#     physics account, at most 40 at once (the M15 runs hold 200 physics + 32 free; limit 240 + 32);
#  3. merge_W14.py -> ../analysis_W14, prefilter_W14.py (KS cut 0.0468, as every leg), initial_vars.csv from
#     ../analysis_W (same var files), ../scmd_noco_pipeline/analyze_scmd_only.py, scmd_5way_compare.py.
# A task is finished when its log has "done" or "already done" (tasks can hang in RUNNING after exiting).
P=/dartfs-hpc/scratch/f007058/scmd_W14_pipeline
S=/dartfs-hpc/scratch/f007058
echo $$ > $P/driver.pid
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate fidanka-env
GEN=$(awk '/gen_W14 array/{print $5}' $S/emul_W14/jobids.txt | tail -1)
[[ "$GEN" =~ ^[0-9]+$ ]] || { echo "no gen_W14 job id in ../emul_W14/jobids.txt"; exit 1; }
echo "=== waiting for isochrones (gen_W14 array $GEN) $(date '+%F %T')"
ngen() { grep -al "^done " $S/emul_W14/logs/gen_${GEN}_*.out 2>/dev/null | wc -l; }
while [ "$(ngen)" -lt 100 ]; do
  squeue -h -j $GEN 2>/dev/null | grep -q . || break
  sleep 600
done
scancel $GEN 2>/dev/null
for d in $S/Isochrones_W14/2*; do
  f=$d/isos/isochrones_example.txt
  [ -s $f ] && [ "$(stat -c %s $f)" -gt 1000000 ] && basename $d
done | sort -n > $P/runs_W14.txt
N=$(wc -l < $P/runs_W14.txt)
echo "isochrones: $N runs (gen tasks done: $(ngen)/100) $(date '+%F %T')"
[ "$N" -gt 0 ] || exit 1

echo "=== sCMD pipeline $(date '+%F %T')"
J=$(sbatch --parsable --array=0-$((N - 1))%40 --export=ALL,RUNLIST=runs_W14.txt $P/run_scmd_W14.batch)
echo "sCMD array $J ($N runs, 40 at once)" | tee -a $P/jobids.txt
ndone() { grep -alE "^done |already done" $P/logs/scmd_${J}_*.out 2>/dev/null | wc -l; }
while [ "$(ndone)" -lt "$N" ]; do
  squeue -h -j $J 2>/dev/null | grep -q . || break
  sleep 900
done
scancel $J 2>/dev/null
echo "sCMD done: $(ndone)/$N tasks, $(ls $S/2d_ks_W14/2d_ks_results_*.csv 2>/dev/null | wc -l) KS result files $(date '+%F %T')"

echo "=== analysis $(date '+%F %T')"
mkdir -p $S/analysis_W14
cp -p $S/analysis_W/initial_vars.csv $S/analysis_W14/
cd $P
python merge_W14.py
python prefilter_W14.py
SCMD_HOME=$S/analysis_W14 SCMD_OUTPUT_DIR=$S/analysis_W14 python $S/scmd_noco_pipeline/analyze_scmd_only.py
python scmd_5way_compare.py
echo "all done $(date '+%F %T')"
