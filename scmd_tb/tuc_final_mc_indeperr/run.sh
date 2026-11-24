#!/bin/bash
# tuc_final_mc_indeperr: the Monte Carlo-on 47 Tuc fit (../tuc_final_mc, with production's refinement override) with
# one change: ast_errors "independent" (each filter's completeness and error drawn from its own radius x magnitude
# bins, as M. Ying's MakeTestCMD; the advisor's model) instead of "nearest" (both filters from one artificial star).
# setup.json = ../tuc_final_mc/setup.json with only ast_errors changed; code/ = a copy of ../tuc_final_mc/code;
# runs.txt, tuc_age_limits_d4mc.json and refine/setup_override.json copied from there / ../tuc_final.
# Steps: 0. wait for ../tuc_final_mc/rerun_refine.sh; 1. grid; 2. KS refinement (E 0.001); 3. bootstrap at this
# fit's own lowest cell; 4. statistics, comparison with ../tuc_final_mc_refined (joint errors) and production.
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export SCMD_PYTHON=$PY HOME=/dartfs-hpc/scratch/f007058 OMP_NUM_THREADS=1
S=/dartfs-hpc/scratch/f007058
B=$S/scmd_tb
K=$B/ks_refine
O=$B/tuc_final_mc_indeperr
V=$B/tuc_final_mc_indeperr_refined
I=$B/tuc_iso_d4v2_mcfix
J=$B/tuc_final_mc
LIM=$O/tuc_age_limits_d4mc.json
export TUC_AGE_LIMITS=$LIM
echo $$ > $O/driver.pid
cd $O/code

echo "=== waiting for the joint-error refinement rerun (../tuc_final_mc/rerun_refine.sh) $(date '+%F %T')"
until grep -q "all done" $J/rerun_refine.log 2>/dev/null || ! kill -0 "$(cat $J/rerun_refine.pid 2>/dev/null)" 2>/dev/null; do sleep 300; done
until ! squeue -h -u f007058 -n tucmc_refine,boot_tucmc | grep -q .; do sleep 120; done
[ -s $O/refine/setup_override.json ] || { echo "no refine/setup_override.json"; exit 1; }
echo "starting $(date '+%F %T'); runs: $(wc -l < $O/runs.txt); $(grep ast_errors $O/setup.json)"
NP=$((240 - $(squeue -u f007058 -h -t R -A physics -o '%C' | awk '{s+=$1} END {print s+0}')))
[ $NP -gt 240 ] && NP=240; [ $NP -lt 100 ] && NP=100
echo "physics workers: $NP (+ 32 free)"

echo "=== grid $(date '+%F %T')"
st() { $PY -m isochronetoolbox.scmd status --out $O --runs $O/runs.txt 2>/dev/null | tail -1; }
G1=$(sbatch --parsable --account=physics --array=0-$((NP - 1)) --job-name=tucind_scmd --output=$O/logs/%A_%a.out scripts/scmd_worker.sbatch $O $O/runs.txt --max-hours 23)
G2=$(sbatch --parsable --account=free --array=0-31 --job-name=tucind_scmd --output=$O/logs/%A_%a.out scripts/scmd_worker.sbatch $O $O/runs.txt --max-hours 23)
echo "grid jobs $G1 $G2 $(date '+%F %T')" | tee -a $O/jobids.txt
until st | grep -q "'claimed': 0, 'todo': 0" || ! squeue -h -u f007058 -n tucind_scmd | grep -q .; do sleep 120; done
sleep 60; scancel -n tucind_scmd -u f007058 2>/dev/null; sleep 20
echo "grid done $(date '+%F %T')"; st
$PY -m isochronetoolbox.scmd collect --out $O 2>&1 | grep -v Warn

echo "=== refinement $(date '+%F %T')"
mkdir -p $O/refine/logs $O/refine/claims
rs() { $PY $K/refine.py $O status 2>/dev/null | tail -1; }
R1=$(sbatch --parsable --export=ALL --account=physics --array=0-$((NP - 1)) --job-name=tucind_refine --output=$O/refine/logs/%A_%a.out $K/refine.sbatch $O $LIM)
R2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --job-name=tucind_refine --output=$O/refine/logs/%A_%a.out $K/refine.sbatch $O $LIM)
echo "refine jobs $R1 $R2" | tee -a $O/jobids.txt
until rs | grep -q "'claimed': 0, 'todo': 0" || ! squeue -h -u f007058 -n tucind_refine | grep -q .; do sleep 120; done
sleep 30; scancel -n tucind_refine -u f007058 2>/dev/null
echo "refinement done $(date '+%F %T')"; rs
$PY $K/refine.py $O merge || exit 1

echo "=== bootstrap $(date '+%F %T')"
mkdir -p $V/logs $V/analysis
cp $O/setup.json $O/runs.txt $V/
$PY -c "import json; s = json.load(open('$V/setup.json')); s.update(json.load(open('$O/refine/setup_override.json'))); json.dump(s, open('$V/setup.json', 'w'), indent=1)" || exit 1
ln -sfn $O/code $V/code
cp $O/combined_2d_ks_results_refined.csv $V/combined_2d_ks_results.csv
CELL=$($PY - <<PYEOF
import json, pandas as pd
lim = {int(k): v for k, v in json.load(open("$LIM")).items()}
g = pd.read_csv("$V/combined_2d_ks_results.csv")
g = g[g.age_gyr < g.simulation_id.map(lim) - 1e-6]
r = g.loc[g.ks_statistic.idxmin()]
print(f"{int(r.simulation_id)}:{r.age_gyr:.4f}:{r.distance_modulus:.4f}:{r.extinction_ebv:.4f}")
PYEOF
)
NT=$((NP + 32))
echo "bootstrap cell $CELL ($NT tasks)"
B1=$(sbatch --parsable --export=ALL --account=physics --array=0-$((NP - 1)) --time=04:00:00 --output=$V/logs/boot_%A_%a.out --job-name=boot_tucind scripts/scmd_bootstrap.sbatch $V --cell $CELL --n 10000 --ntasks $NT --task-offset 0 --fresh-reference)
B2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --time=04:00:00 --output=$V/logs/boot_%A_%a.out --job-name=boot_tucind scripts/scmd_bootstrap.sbatch $V --cell $CELL --n 10000 --ntasks $NT --task-offset $NP --fresh-reference)
echo "bootstrap jobs $B1 $B2" | tee -a $V/jobids.txt $O/jobids.txt
until [ "$(ls $V/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge $NT ]; do
    squeue -h -u f007058 -n boot_tucind | grep -q . || break
    sleep 60
done
sleep 30; scancel -n boot_tucind -u f007058 2>/dev/null
echo "bootstrap done $(date '+%F %T'): $(ls $V/bootstrap/part_*.csv 2>/dev/null | wc -l) parts"

echo "=== statistics $(date '+%F %T')"
$PY -m isochronetoolbox.scmd bootstrap-collect --out $V 2>&1 | grep -v Warn
echo "--- whole grid: mu 13.15-13.55, E 0.000-0.050"
$PY $B/tucf_statistics.py $V 2>&1 | grep -v Warn
echo "--- distance prior (D17): mu 13.15-13.45, E 0-0.05"
MU_MIN=13.15 MU_MAX=13.45 EVI_MIN=0 EVI_MAX=0.05 POST_SUFFIX=_D17 $PY $B/tucf_statistics.py $V 2>&1 | grep -v Warn
echo "--- box converted: mu 13.15-13.45, E 0-0.039"
MU_MIN=13.15 MU_MAX=13.45 EVI_MIN=0 EVI_MAX=0.039 POST_SUFFIX=_D17conv $PY $B/tucf_statistics.py $V 2>&1 | grep -v Warn
$PY $B/best3_and_null.py $V "DSEP4, Monte Carlo on, independent AST errors" 2>&1 | grep -v Warn
$PY $B/refined_plots.py $V $LIM "47 Tuc, Monte Carlo on, independent AST errors" 2>&1 | grep -v Warn
echo "--- age posteriors: production (Monte Carlo off), joint errors (Monte Carlo on), independent errors (Monte Carlo on)"
for p in strict strict_D17 shifted shifted_D17 scaled scaled_D17; do
  for g in $B/tuc_final_refined $B/tuc_final_mc_refined $V; do
    f=$g/posterior_$p/posterior_summary.csv
    [ -s $f ] && awk -F, -v p=$p -v g=$(basename $g) '$1=="age" {printf "%-22s %-24s mean %.2f median %.2f sd %.2f  16-84%% %.2f-%.2f  runs %s\n", p, g, $2, $3, $7, $5, $6, $9}' $f
  done
done | tee $O/compare_joint_vs_independent.txt
echo "all done $(date '+%F %T')"
