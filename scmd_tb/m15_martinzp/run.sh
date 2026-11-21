#!/bin/bash
# m15_martinzp: the production M15 pipeline (m15_final -> m15_final_refined) on the same 500 runs with one
# change: M. Ying's colour zero-point correction (M15_ZPCR.ipynb, ast.martin_zero_point) in place of ours.
# Everything else, the statistics included, is production: ridge cut, nearest-AST errors, raw radii, PDMF
# -1.00, binaries 0.017, stage 2 + 3 grid, KS refinement (E step 0.001, E 0.000-0.120), 10,000-realisation
# bootstrap (fresh references, resampled ASTs) at the lowest-KS cell, strict / shifted / scaled posteriors.
#  1. grid: physics 240 + free 32 workers
#  2. collect; KS refinement (ks_refine/refine.py, same workers); merge
#  3. m15_martinzp_refined/: setup with the E override, code link, refined rows; bootstrap (272 tasks)
#  4. statistics (all E; E >= 0.05; E >= 0.08; the Bhardwaj/Dutra & Bica box, as is and x0.78), figures
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export SCMD_PYTHON=$PY HOME=/dartfs-hpc/scratch/f007058 OMP_NUM_THREADS=1
B=/dartfs-hpc/scratch/f007058/scmd_tb
K=$B/ks_refine
O=$B/m15_martinzp
V=$B/m15_martinzp_refined
LIM=/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json
export TUC_AGE_LIMITS=$LIM
echo $$ > $O/driver.pid
cd $O/code

echo "=== grid $(date '+%F %T')"
st() { $PY -m isochronetoolbox.scmd status --out $O --runs $O/runs.txt 2>/dev/null | tail -1; }
G1=$(sbatch --parsable --account=physics --array=0-239 --job-name=m15mz_scmd --output=$O/logs/%A_%a.out scripts/scmd_worker.sbatch $O $O/runs.txt --max-hours 23)
G2=$(sbatch --parsable --account=free --array=0-31 --job-name=m15mz_scmd --output=$O/logs/%A_%a.out scripts/scmd_worker.sbatch $O $O/runs.txt --max-hours 23)
echo "grid jobs $G1 $G2 $(date '+%F %T')" | tee -a $O/jobids.txt
until st | grep -q "'claimed': 0, 'todo': 0" || ! squeue -h -u f007058 -n m15mz_scmd | grep -q .; do sleep 120; done
sleep 60; scancel -n m15mz_scmd -u f007058 2>/dev/null; sleep 20
echo "grid done $(date '+%F %T')"; st
$PY -m isochronetoolbox.scmd collect --out $O 2>&1 | grep -v Warn

echo "=== refinement $(date '+%F %T')"
mkdir -p $O/refine/logs $O/refine/claims
rs() { $PY $K/refine.py $O status 2>/dev/null | tail -1; }
R1=$(sbatch --parsable --export=ALL --account=physics --array=0-239 --job-name=m15mz_refine --output=$O/refine/logs/%A_%a.out $K/refine.sbatch $O $LIM)
R2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --job-name=m15mz_refine --output=$O/refine/logs/%A_%a.out $K/refine.sbatch $O $LIM)
echo "refine jobs $R1 $R2" | tee -a $O/jobids.txt
until rs | grep -q "'claimed': 0, 'todo': 0" || ! squeue -h -u f007058 -n m15mz_refine | grep -q .; do sleep 120; done
sleep 30; scancel -n m15mz_refine -u f007058 2>/dev/null
echo "refinement done $(date '+%F %T')"; rs
$PY $K/refine.py $O merge || exit 1

echo "=== bootstrap $(date '+%F %T')"
mkdir -p $V/logs $V/analysis
cp $O/setup.json $O/runs.txt $V/
$PY -c "import json; s = json.load(open('$V/setup.json')); s.update(json.load(open('$O/refine/setup_override.json'))); json.dump(s, open('$V/setup.json', 'w'), indent=1)"
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
echo "bootstrap cell $CELL"
B1=$(sbatch --parsable --export=ALL --account=physics --array=0-239 --time=03:00:00 --output=$V/logs/boot_%A_%a.out --job-name=boot_m15mz scripts/scmd_bootstrap.sbatch $V --cell $CELL --n 10000 --ntasks 272 --task-offset 0 --fresh-reference)
B2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --time=03:00:00 --output=$V/logs/boot_%A_%a.out --job-name=boot_m15mz scripts/scmd_bootstrap.sbatch $V --cell $CELL --n 10000 --ntasks 272 --task-offset 240 --fresh-reference)
echo "bootstrap jobs $B1 $B2" | tee -a $V/jobids.txt $O/jobids.txt
until [ "$(ls $V/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 272 ]; do
    squeue -h -u f007058 -n boot_m15mz | grep -q . || break
    sleep 60
done
sleep 30; scancel -n boot_m15mz -u f007058 2>/dev/null
echo "bootstrap done $(date '+%F %T'): $(ls $V/bootstrap/part_*.csv 2>/dev/null | wc -l) parts"

echo "=== statistics $(date '+%F %T')"
$PY -m isochronetoolbox.scmd bootstrap-collect --out $V 2>&1 | grep -v Warn
echo "--- all E (0.000-0.120)"
$PY $B/m15f_statistics.py $V 2>&1 | grep -v Warn
echo "--- E >= 0.05"
EVI_MIN=0.05 POST_SUFFIX=_E05 $PY $B/m15f_statistics.py $V 2>&1 | grep -v Warn
echo "--- E >= 0.08"
EVI_MIN=0.08 POST_SUFFIX=_Elit $PY $B/m15f_statistics.py $V 2>&1 | grep -v Warn
echo "--- box mu 15.3-15.6, E 0.08-0.15"
MU_MIN=15.3 MU_MAX=15.6 EVI_MIN=0.08 EVI_MAX=0.15 POST_SUFFIX=_DB $PY $B/m15f_statistics.py $V 2>&1 | grep -v Warn
echo "--- box mu 15.3-15.6, E 0.062-0.117 (E(V-I) x 0.78)"
MU_MIN=15.3 MU_MAX=15.6 EVI_MIN=0.062 EVI_MAX=0.117 POST_SUFFIX=_DBconv $PY $B/m15f_statistics.py $V 2>&1 | grep -v Warn
$PY $B/m15_best3_and_null.py $V DSEP4 2>&1 | grep -v Warn
$PY $B/refined_plots.py $V $LIM "M15, Martin zero point" 2>&1 | grep -v Warn
echo "all done $(date '+%F %T')"
