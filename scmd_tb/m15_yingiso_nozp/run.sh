#!/bin/bash
# m15_yingiso_nozp: recreation of Ying et al. (2025)'s M15 fit with his published isochrones (Zenodo 13998637,
# converted to ../../M15/iso_ying2025 by convert_ying_isochrones.py; 41 ages, 8-16 Gyr in 0.2 Gyr steps). Everything
# else is ../m15_yingbox_nozp (his logic where we have it): no spatial zero-point map, his E box 0.08-0.15 applied
# directly, mu 15.30-15.60; our sCMDs (ridge cut, nearest-AST errors, PDMF -1.00, binaries 0.017) and our
# 95th-percentile 2D KS, grid + KS refinement, 10,000-realisation bootstrap (fresh references) at the lowest cell,
# strict / shifted / scaled posteriors. The one change from ../m15_yingbox_nozp is the isochrones.
#  1. grid: physics 240 + free 32 workers; collect
#  2. KS refinement (ks_refine/refine.py, same workers); merge
#  3. ../m15_yingiso_nozp_refined/: setup with the E override, code link, refined rows; bootstrap (272 tasks)
#  4. statistics (m15yi_statistics.py: bins to 16 Gyr and E 0.155), figures, comparison with our isochrones
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export SCMD_PYTHON=$PY HOME=/dartfs-hpc/scratch/f007058 OMP_NUM_THREADS=1
B=/dartfs-hpc/scratch/f007058/scmd_tb
K=$B/ks_refine
O=$B/m15_yingiso_nozp
V=$B/m15_yingiso_nozp_refined
LIM=/dartfs-hpc/scratch/f007058/M15/iso_ying2025/age_limits_ying.json
export TUC_AGE_LIMITS=$LIM
echo $$ > $O/driver.pid
cd $O/code

echo "=== grid $(date '+%F %T')"
st() { $PY -m isochronetoolbox.scmd status --out $O --runs $O/runs.txt 2>/dev/null | tail -1; }
G1=$(sbatch --parsable --account=physics --array=0-239 --job-name=m15yi_scmd --output=$O/logs/%A_%a.out scripts/scmd_worker.sbatch $O $O/runs.txt --max-hours 23)
G2=$(sbatch --parsable --account=free --array=0-31 --job-name=m15yi_scmd --output=$O/logs/%A_%a.out scripts/scmd_worker.sbatch $O $O/runs.txt --max-hours 23)
echo "grid jobs $G1 $G2 $(date '+%F %T')" | tee -a $O/jobids.txt
until st | grep -q "'claimed': 0, 'todo': 0" || ! squeue -h -u f007058 -n m15yi_scmd | grep -q .; do sleep 120; done
sleep 60; scancel -n m15yi_scmd -u f007058 2>/dev/null; sleep 20
echo "grid done $(date '+%F %T')"; st
$PY -m isochronetoolbox.scmd collect --out $O 2>&1 | grep -v Warn

echo "=== refinement $(date '+%F %T')"
mkdir -p $O/refine/logs $O/refine/claims
rs() { $PY $K/refine.py $O status 2>/dev/null | tail -1; }
R1=$(sbatch --parsable --export=ALL --account=physics --array=0-239 --job-name=m15yi_refine --output=$O/refine/logs/%A_%a.out $K/refine.sbatch $O $LIM)
R2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --job-name=m15yi_refine --output=$O/refine/logs/%A_%a.out $K/refine.sbatch $O $LIM)
echo "refine jobs $R1 $R2" | tee -a $O/jobids.txt
until rs | grep -q "'claimed': 0, 'todo': 0" || ! squeue -h -u f007058 -n m15yi_refine | grep -q .; do sleep 120; done
sleep 30; scancel -n m15yi_refine -u f007058 2>/dev/null
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
B1=$(sbatch --parsable --export=ALL --account=physics --array=0-239 --time=03:00:00 --output=$V/logs/boot_%A_%a.out --job-name=boot_m15yi scripts/scmd_bootstrap.sbatch $V --cell $CELL --n 10000 --ntasks 272 --task-offset 0 --fresh-reference)
B2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --time=03:00:00 --output=$V/logs/boot_%A_%a.out --job-name=boot_m15yi scripts/scmd_bootstrap.sbatch $V --cell $CELL --n 10000 --ntasks 272 --task-offset 240 --fresh-reference)
echo "bootstrap jobs $B1 $B2" | tee -a $V/jobids.txt $O/jobids.txt
until [ "$(ls $V/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 272 ]; do
    squeue -h -u f007058 -n boot_m15yi | grep -q . || break
    sleep 60
done
sleep 30; scancel -n boot_m15yi -u f007058 2>/dev/null
echo "bootstrap done $(date '+%F %T'): $(ls $V/bootstrap/part_*.csv 2>/dev/null | wc -l) parts"

echo "=== statistics $(date '+%F %T')"
$PY -m isochronetoolbox.scmd bootstrap-collect --out $V 2>&1 | grep -v Warn
echo "--- Ying's box: mu 15.30-15.60, E 0.080-0.150 (the whole grid), Ying's isochrones"
$PY $O/m15yi_statistics.py $V 2>&1 | grep -v Warn
$PY $B/m15_best3_and_null.py $V "Ying iso" 2>&1 | grep -v Warn
$PY $B/refined_plots.py $V $LIM "M15, Ying's isochrones, Ying's E box, no zero-point map" 2>&1 | grep -v Warn
$PY $O/compare_ks_vs_age.py 2>&1 | grep -v Warn
echo "all done $(date '+%F %T')"
