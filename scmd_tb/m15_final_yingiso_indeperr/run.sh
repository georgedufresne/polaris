#!/bin/bash
# m15_final_yingiso_indeperr: production (../m15_final: zero-point map, E 0-0.12, ridge cut, 95th-percentile KS, grid + KS
# refinement, bootstrap, strict / shifted / scaled posteriors) with Ying's published isochrones (Zenodo 13998637; his
# 8-16 Gyr grid in 0.2) and the advisor's independent per-filter artificial-star errors (ast_errors "independent").
# Pairs: ../m15_final_indeperr (same, our isochrones) isolates the isochrones; ../m15_final is production.
#  1. grid: physics 240 + free 32 workers; collect
#  2. KS refinement (ks_refine/refine.py, same workers); merge
#  3. ../m15_final_yingiso_indeperr_refined/: setup with the E override, code link, refined rows; bootstrap (272 tasks)
#  4. statistics (m15yi_statistics.py: bins to 16 Gyr), figures, comparison with the other runs
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export SCMD_PYTHON=$PY HOME=/dartfs-hpc/scratch/f007058 OMP_NUM_THREADS=1
B=/dartfs-hpc/scratch/f007058/scmd_tb
K=$B/ks_refine
O=$B/m15_final_yingiso_indeperr
V=$B/m15_final_yingiso_indeperr_refined
LIM=/dartfs-hpc/scratch/f007058/M15/iso_ying2025/age_limits_ying.json
export TUC_AGE_LIMITS=$LIM
echo $$ > $O/driver.pid
cd $O/code

echo "=== grid $(date '+%F %T')"
st() { $PY -m isochronetoolbox.scmd status --out $O --runs $O/runs.txt 2>/dev/null | tail -1; }
G1=$(sbatch --parsable --account=physics --array=0-239 --job-name=m15yie_scmd --output=$O/logs/%A_%a.out scripts/scmd_worker.sbatch $O $O/runs.txt --max-hours 23)
G2=$(sbatch --parsable --account=free --array=0-31 --job-name=m15yie_scmd --output=$O/logs/%A_%a.out scripts/scmd_worker.sbatch $O $O/runs.txt --max-hours 23)
echo "grid jobs $G1 $G2 $(date '+%F %T')" | tee -a $O/jobids.txt
until st | grep -q "'claimed': 0, 'todo': 0" || ! squeue -h -u f007058 -n m15yie_scmd | grep -q .; do sleep 120; done
sleep 60; scancel -n m15yie_scmd -u f007058 2>/dev/null; sleep 20
echo "grid done $(date '+%F %T')"; st
$PY -m isochronetoolbox.scmd collect --out $O 2>&1 | grep -v Warn

echo "=== refinement $(date '+%F %T')"
mkdir -p $O/refine/logs $O/refine/claims
rs() { $PY $K/refine.py $O status 2>/dev/null | tail -1; }
R1=$(sbatch --parsable --export=ALL --account=physics --array=0-239 --job-name=m15yie_refine --output=$O/refine/logs/%A_%a.out $K/refine.sbatch $O $LIM)
R2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --job-name=m15yie_refine --output=$O/refine/logs/%A_%a.out $K/refine.sbatch $O $LIM)
echo "refine jobs $R1 $R2" | tee -a $O/jobids.txt
until rs | grep -q "'claimed': 0, 'todo': 0" || ! squeue -h -u f007058 -n m15yie_refine | grep -q .; do sleep 120; done
sleep 30; scancel -n m15yie_refine -u f007058 2>/dev/null
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
B1=$(sbatch --parsable --export=ALL --account=physics --array=0-239 --time=03:00:00 --output=$V/logs/boot_%A_%a.out --job-name=boot_m15yie scripts/scmd_bootstrap.sbatch $V --cell $CELL --n 10000 --ntasks 272 --task-offset 0 --fresh-reference)
B2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --time=03:00:00 --output=$V/logs/boot_%A_%a.out --job-name=boot_m15yie scripts/scmd_bootstrap.sbatch $V --cell $CELL --n 10000 --ntasks 272 --task-offset 240 --fresh-reference)
echo "bootstrap jobs $B1 $B2" | tee -a $V/jobids.txt $O/jobids.txt
until [ "$(ls $V/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 272 ]; do
    squeue -h -u f007058 -n boot_m15yie | grep -q . || break
    sleep 60
done
sleep 30; scancel -n boot_m15yie -u f007058 2>/dev/null
echo "bootstrap done $(date '+%F %T'): $(ls $V/bootstrap/part_*.csv 2>/dev/null | wc -l) parts"

echo "=== statistics $(date '+%F %T')"
$PY -m isochronetoolbox.scmd bootstrap-collect --out $V 2>&1 | grep -v Warn
echo "--- production setup with Ying's isochrones and independent per-filter errors"
$PY $O/m15yi_statistics.py $V 2>&1 | grep -v Warn
$PY $B/m15_best3_and_null.py $V "production + Ying iso + independent errors" 2>&1 | grep -v Warn
$PY $B/refined_plots.py $V $LIM "M15, production with Ying's isochrones and independent per-filter errors" 2>&1 | grep -v Warn
$PY $O/compare_runs.py 2>&1 | grep -v Warn | tee $O/logs/compare_runs.out
echo "all done $(date '+%F %T')"
