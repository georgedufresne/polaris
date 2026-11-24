#!/bin/bash
# Rerun of run.sh steps 2-4 (refinement -> bootstrap -> statistics) with production's refinement override
# (refine/setup_override.json = ../tuc_final/refine/setup_override.json: E 0.000-0.050 in 0.001). The first pass
# (2026-11-24) had no override file, so the descent ran on the coarse 0.005 E grid; its outputs are kept as
# refine_coarseE/, combined_2d_ks_results_refined_coarseE.csv, ../tuc_final_mc_refined_coarseE/.
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export SCMD_PYTHON=$PY HOME=/dartfs-hpc/scratch/f007058 OMP_NUM_THREADS=1
S=/dartfs-hpc/scratch/f007058
B=$S/scmd_tb
K=$B/ks_refine
O=$B/tuc_final_mc
V=$B/tuc_final_mc_refined
LIM=$O/tuc_age_limits_d4mc.json
export TUC_AGE_LIMITS=$LIM
echo $$ > $O/rerun_refine.pid
cd $O/code
[ -s $O/refine/setup_override.json ] || { echo "no refine/setup_override.json"; exit 1; }
NP=$((240 - $(squeue -u f007058 -h -t R -A physics -o '%C' | awk '{s+=$1} END {print s+0}')))
[ $NP -gt 240 ] && NP=240; [ $NP -lt 100 ] && NP=100
echo "physics workers: $NP (+ 32 free); override $(cat $O/refine/setup_override.json)"
echo "=== refinement $(date '+%F %T')"
mkdir -p $O/refine/logs $O/refine/claims
rs() { $PY $K/refine.py $O status 2>/dev/null | tail -1; }
R1=$(sbatch --parsable --export=ALL --account=physics --array=0-$((NP - 1)) --job-name=tucmc_refine --output=$O/refine/logs/%A_%a.out $K/refine.sbatch $O $LIM)
R2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --job-name=tucmc_refine --output=$O/refine/logs/%A_%a.out $K/refine.sbatch $O $LIM)
echo "refine jobs $R1 $R2" | tee -a $O/jobids.txt
until rs | grep -q "'claimed': 0, 'todo': 0" || ! squeue -h -u f007058 -n tucmc_refine | grep -q .; do sleep 120; done
sleep 30; scancel -n tucmc_refine -u f007058 2>/dev/null
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
B1=$(sbatch --parsable --export=ALL --account=physics --array=0-$((NP - 1)) --time=04:00:00 --output=$V/logs/boot_%A_%a.out --job-name=boot_tucmc scripts/scmd_bootstrap.sbatch $V --cell $CELL --n 10000 --ntasks $NT --task-offset 0 --fresh-reference)
B2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --time=04:00:00 --output=$V/logs/boot_%A_%a.out --job-name=boot_tucmc scripts/scmd_bootstrap.sbatch $V --cell $CELL --n 10000 --ntasks $NT --task-offset $NP --fresh-reference)
echo "bootstrap jobs $B1 $B2" | tee -a $V/jobids.txt $O/jobids.txt
until [ "$(ls $V/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge $NT ]; do
    squeue -h -u f007058 -n boot_tucmc | grep -q . || break
    sleep 60
done
sleep 30; scancel -n boot_tucmc -u f007058 2>/dev/null
echo "bootstrap done $(date '+%F %T'): $(ls $V/bootstrap/part_*.csv 2>/dev/null | wc -l) parts"

echo "=== statistics $(date '+%F %T')"
$PY -m isochronetoolbox.scmd bootstrap-collect --out $V 2>&1 | grep -v Warn
echo "--- whole grid: mu 13.15-13.55, E 0.000-0.050"
$PY $B/tucf_statistics.py $V 2>&1 | grep -v Warn
echo "--- distance prior (D17): mu 13.15-13.45, E 0-0.05"
MU_MIN=13.15 MU_MAX=13.45 EVI_MIN=0 EVI_MAX=0.05 POST_SUFFIX=_D17 $PY $B/tucf_statistics.py $V 2>&1 | grep -v Warn
echo "--- box converted: mu 13.15-13.45, E 0-0.039"
MU_MIN=13.15 MU_MAX=13.45 EVI_MIN=0 EVI_MAX=0.039 POST_SUFFIX=_D17conv $PY $B/tucf_statistics.py $V 2>&1 | grep -v Warn
$PY $B/best3_and_null.py $V "DSEP4, Monte Carlo on" 2>&1 | grep -v Warn
$PY $B/refined_plots.py $V $LIM "47 Tuc, production on the Monte Carlo-on DSEP4 isochrones" 2>&1 | grep -v Warn
echo "--- age posteriors: production (Monte Carlo off) vs this (Monte Carlo on)"
for p in strict strict_D17 shifted shifted_D17 scaled scaled_D17; do
  for g in $B/tuc_final_refined $V; do
    f=$g/posterior_$p/posterior_summary.csv
    [ -s $f ] && awk -F, -v p=$p -v g=$(basename $g) '$1=="age" {printf "%-22s %-24s mean %.2f median %.2f sd %.2f  16-84%% %.2f-%.2f  runs %s\n", p, g, $2, $3, $7, $5, $6, $9}' $f
  done
done | tee $O/compare_production.txt; for p in strict strict_D17 shifted_D17 scaled_D17; do f=$B/tuc_final_mc_refined_coarseE/posterior_$p/posterior_summary.csv; [ -s $f ] && awk -F, -v p=$p '$1=="age" {printf "%-22s %-24s mean %.2f median %.2f sd %.2f  runs %s\n", p, "coarse-E first pass", $2, $3, $7, $9}' $f; done | tee -a $O/compare_production.txt
echo "all done $(date '+%F %T')"
