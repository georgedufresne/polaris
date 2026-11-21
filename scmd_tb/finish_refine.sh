#!/bin/bash
# KS refinement of m15_final and tuc_final (ks_refine/refine.py), then bootstraps and posteriors on the
# refined grids. Stage 2 (KDE JS, keep 12) passed cells to stage 3 that are far from the 2D-KS minimum
# (m15_final run 10392, 13.6 Gyr: kept KS 0.050, true minimum 0.0150 at a cell ranked 47th by stage 2).
#  1. finish the 47 Tuc grid (as finish_tuc_final.sh, stopped before its bootstrap), collect
#  2. refinement workers: physics 240 + free 32, M15 units first (its 8 chain runs first), then 47 Tuc
#  3. validate against the chains' exhaustive scoring, merge -> <grid>/combined_2d_ks_results_refined.csv
#  4. <grid>_refined/: setup, code link, refined rows; bootstraps there (10,000, fresh references, resampled
#     ASTs, 136 tasks each, run together) at the lowest-KS cell below the age limits
#  5. statistics (M15: all E, E >= 0.05 = the old grid floor, E >= 0.08 = literature E(B-V) 0.10 +/- 0.02),
#     best-3 figures
# The descent runs on E steps of 0.001 (<grid>/refine/setup_override.json; the KS valley in age-E is
# narrower than the grid's 0.005), and for M15 down to E 0.000 (E(B-V) prior lowered, 2026-11-05). The
# _refined/setup.json carries the override so bootstraps and replays can regenerate those cells.
# E is E(F606W-F814W) = 1.02 E(B-V) (ACS 2.88 - 1.86), mu the apparent F606W modulus.
# The original grids, their combined_2d_ks_results.csv and m15_final's bootstrap are not modified.
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export SCMD_PYTHON=$PY HOME=/dartfs-hpc/scratch/f007058 OMP_NUM_THREADS=1
B=/dartfs-hpc/scratch/f007058/scmd_tb
K=$B/ks_refine
M=$B/m15_final
T=$B/tuc_final
LM=/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json
LT=$T/tuc_age_limits_d4.json
mkdir -p $K/logs
echo $$ > $K/driver.pid

echo "=== 1. 47 Tuc grid $(date '+%F %T')"
cd $T/code
st() { $PY -m isochronetoolbox.scmd status --out $T --runs $T/runs.txt 2>/dev/null | grep "'done'"; }
until st | grep -q "'claimed': 0, 'todo': 0" || ! squeue -h -u f007058 -n tucf_scmd | grep -q .; do sleep 120; done
sleep 60; scancel -n tucf_scmd -u f007058 2>/dev/null; sleep 20
echo "grid done $(date '+%F %T')"; st
$PY -m isochronetoolbox.scmd collect --out $T 2>&1 | grep -v Warn

echo "=== 2. refinement $(date '+%F %T')"
cd $K
export REFINE_FIRST_1=$K/first_m15.txt REFINE_FIRST_2=$K/first_tuc.txt
R1=$(sbatch --parsable --export=ALL --account=physics --array=0-239 --output=$K/logs/%A_%a.out refine.sbatch $M $LM $T $LT)
R2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --output=$K/logs/%A_%a.out refine.sbatch $M $LM $T $LT)
echo "refine jobs $R1 $R2" | tee -a $M/jobids.txt $T/jobids.txt
rs() { $PY $K/refine.py $1 status 2>/dev/null | tail -1; }
until { rs $M | grep -q "'claimed': 0, 'todo': 0" && rs $T | grep -q "'claimed': 0, 'todo': 0"; } \
      || ! squeue -h -u f007058 -n ks_refine | grep -q .; do sleep 120; done
sleep 30; scancel -n ks_refine -u f007058 2>/dev/null
echo "refinement done $(date '+%F %T')"; rs $M; rs $T

echo "=== 3. validation and merge"
$PY $K/validate.py m15
$PY $K/validate.py tuc
$PY $K/refine.py $M merge || exit 1
$PY $K/refine.py $T merge || exit 1

echo "=== 4. bootstraps $(date '+%F %T')"
declare -A LIM=([m15]=$LM [tuc]=$LT)
for X in m15 tuc; do
    G=$B/${X}_final; V=$B/${X}_final_refined
    mkdir -p $V/logs $V/analysis
    cp $G/setup.json $G/runs.txt $V/
    if [ -f $G/refine/setup_override.json ]; then
        $PY -c "import json; s = json.load(open('$V/setup.json')); s.update(json.load(open('$G/refine/setup_override.json'))); json.dump(s, open('$V/setup.json', 'w'), indent=1)"
        echo "$X setup override: $(cat $G/refine/setup_override.json)"
    fi
    ln -sfn $G/code $V/code
    cp $G/combined_2d_ks_results_refined.csv $V/combined_2d_ks_results.csv
    CELL=$($PY - <<PYEOF
import json, pandas as pd
lim = {int(k): v for k, v in json.load(open("${LIM[$X]}")).items()}
g = pd.read_csv("$V/combined_2d_ks_results.csv")
g = g[g.age_gyr < g.simulation_id.map(lim) - 1e-6]
r = g.loc[g.ks_statistic.idxmin()]
print(f"{int(r.simulation_id)}:{r.age_gyr:.4f}:{r.distance_modulus:.4f}:{r.extinction_ebv:.4f}")
PYEOF
)
    echo "$X bootstrap cell $CELL"
    cd $G/code
    B1=$(sbatch --parsable --export=ALL --account=physics --array=0-119 --time=02:00:00 --output=$V/logs/boot_%A_%a.out --job-name=boot_ref scripts/scmd_bootstrap.sbatch $V --cell $CELL --n 10000 --ntasks 136 --task-offset 0 --fresh-reference)
    B2=$(sbatch --parsable --export=ALL --account=free --array=0-15 --time=02:00:00 --output=$V/logs/boot_%A_%a.out --job-name=boot_ref scripts/scmd_bootstrap.sbatch $V --cell $CELL --n 10000 --ntasks 136 --task-offset 120 --fresh-reference)
    echo "$X bootstrap jobs $B1 $B2" | tee -a $V/jobids.txt
done
until [ "$(ls $B/m15_final_refined/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 136 ] \
      && [ "$(ls $B/tuc_final_refined/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 136 ]; do
    squeue -h -u f007058 -n boot_ref | grep -q . || break
    sleep 60
done
sleep 30; scancel -n boot_ref -u f007058 2>/dev/null
echo "bootstraps done $(date '+%F %T')"

echo "=== 5. statistics"
VM=$B/m15_final_refined; VT=$B/tuc_final_refined
cd $M/code; $PY -m isochronetoolbox.scmd bootstrap-collect --out $VM 2>&1 | grep -v Warn
cd $T/code; $PY -m isochronetoolbox.scmd bootstrap-collect --out $VT 2>&1 | grep -v Warn
echo "--- M15, all E (0.000-0.120)"
$PY $B/m15f_statistics.py $VM 2>&1 | grep -v Warn
echo "--- M15, E >= 0.05 (the original grid floor)"
EVI_MIN=0.05 POST_SUFFIX=_E05 $PY $B/m15f_statistics.py $VM 2>&1 | grep -v Warn
echo "--- M15, E >= 0.08 (literature E(B-V) 0.10 +/- 0.02)"
EVI_MIN=0.08 POST_SUFFIX=_Elit $PY $B/m15f_statistics.py $VM 2>&1 | grep -v Warn
TUC_AGE_LIMITS=$LM $PY $B/m15_best3_and_null.py $VM DSEP4 2>&1 | grep -v Warn
echo "--- 47 Tuc"
$PY $B/tucf_statistics.py $VT 2>&1 | grep -v Warn
TUC_AGE_LIMITS=$LT $PY $B/best3_and_null.py $VT DSEP4 2>&1 | grep -v Warn
echo "all done $(date '+%F %T')"
