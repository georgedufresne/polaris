#!/bin/bash
# Paired comparison on this experiment's 500 runs (runs.txt = the first 500 of ../tuc_final_mc/runs.txt, 10000-10500).
# 1. After ../tuc_final_mc/rerun_refine.sh: the joint-error (Monte Carlo on, E refined in 0.001) rows of those 500
#    runs -> joint500/, statistics with the joint fit's own bootstrap null (whole, _D17, _D17conv).
# 2. After this experiment's run.sh: compare_500.txt = age posteriors of production (996, Monte Carlo off), joint (996),
#    joint (500) and independent (500).
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export HOME=/dartfs-hpc/scratch/f007058 OMP_NUM_THREADS=1
B=/dartfs-hpc/scratch/f007058/scmd_tb
J=$B/tuc_final_mc; R=$B/tuc_final_mc_refined; O=$B/tuc_final_mc_indeperr; V=$B/tuc_final_mc_indeperr_refined; G=$O/joint500
echo $$ > $O/joint500.pid
until grep -q "all done" $J/rerun_refine.log 2>/dev/null || ! kill -0 "$(cat $J/rerun_refine.pid 2>/dev/null)" 2>/dev/null; do sleep 300; done
grep -q "all done" $J/rerun_refine.log || { echo "rerun_refine.sh did not finish"; exit 1; }
mkdir -p $G/bootstrap
$PY - <<PYEOF || exit 1
import pandas as pd, numpy as np
runs = set(int(x) for x in open("$O/runs.txt").read().split())
r = pd.read_csv("$R/combined_2d_ks_results.csv")
ne = np.round(r.extinction_ebv, 4).nunique()
assert ne > 40, f"refined rows have only {ne} E values: the 0.001 override is missing"
s = r[r.simulation_id.isin(runs)]
s.to_csv("$G/combined_2d_ks_results.csv", index=False)
print(f"joint rows: {len(r):,} (996 runs, {ne} E values) -> {len(s):,} rows of {s.simulation_id.nunique()} runs")
PYEOF
cp $R/bootstrap_2d_ks_results.csv $R/setup.json $O/runs.txt $G/; cp $R/bootstrap/cell.json $G/bootstrap/
echo "--- joint 500: whole grid"; $PY $B/tucf_statistics.py $G 2>&1 | grep -v Warn
echo "--- joint 500: D17"; MU_MIN=13.15 MU_MAX=13.45 EVI_MIN=0 EVI_MAX=0.05 POST_SUFFIX=_D17 $PY $B/tucf_statistics.py $G 2>&1 | grep -v Warn
echo "--- joint 500: D17conv"; MU_MIN=13.15 MU_MAX=13.45 EVI_MIN=0 EVI_MAX=0.039 POST_SUFFIX=_D17conv $PY $B/tucf_statistics.py $G 2>&1 | grep -v Warn
echo "joint500 statistics done $(date '+%F %T')"
until grep -q "all done" $O/run.log 2>/dev/null || ! kill -0 "$(cat $O/driver.pid 2>/dev/null)" 2>/dev/null; do sleep 600; done
for p in strict strict_D17 shifted shifted_D17 scaled scaled_D17; do
  for g in $B/tuc_final_refined:production_996_MCoff $R:joint_996 $G:joint_500 $V:independent_500; do
    f=${g%%:*}/posterior_$p/posterior_summary.csv
    [ -s $f ] && awk -F, -v p=$p -v g=${g##*:} '$1=="age" {printf "%-12s %-22s mean %.2f median %.2f sd %.2f  16-84%% %.2f-%.2f  runs %s\n", p, g, $2, $3, $7, $5, $6, $9}' $f
  done
done | tee $O/compare_500.txt
echo "all done $(date '+%F %T')"
