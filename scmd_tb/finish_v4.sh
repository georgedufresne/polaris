#!/bin/bash
# tuc_grid_v4 (v2 without the zero-point correction and the PSF scatter): collect,
# bootstrap (fresh references, resampled ASTs) at the best cell within track coverage,
# statistics, best-3 figure.
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export SCMD_PYTHON=$PY HOME=/dartfs-hpc/scratch/f007058 OMP_NUM_THREADS=1
B=/dartfs-hpc/scratch/f007058/scmd_tb
O=$B/tuc_grid_v4
cd $O/code
until [ "$($PY -m isochronetoolbox.scmd status --out $O --runs $O/runs.txt 2>/dev/null | grep -o "'done': [0-9]*" | grep -o "[0-9]*$")" = "679" ]; do sleep 60; done
sleep 60; scancel -n tuc_v4 -u f007058 2>/dev/null; sleep 20
echo "grid done $(date '+%F %T')"
$PY -m isochronetoolbox.scmd collect --out $O 2>&1 | grep -v Warn
CELL=$($PY - <<PYEOF
import json, pandas as pd
lim = {int(k): v for k, v in json.load(open("$B/tuc_age_limits.json")).items()}
g = pd.read_csv("$O/combined_2d_ks_results.csv")
g = g[g.age_gyr < g.simulation_id.map(lim) - 1e-6]
r = g.loc[g.ks_statistic.idxmin()]
print(f"{int(r.simulation_id)}:{r.age_gyr:.4f}:{r.distance_modulus:.4f}:{r.extinction_ebv:.4f}")
PYEOF
)
echo "bootstrap cell $CELL"
B1=$(sbatch --parsable --export=ALL --account=physics --array=0-239 --time=03:00:00 --output=$O/logs/boot_%A_%a.out --job-name=boot_v4 scripts/scmd_bootstrap.sbatch $O --cell $CELL --n 10000 --ntasks 272 --task-offset 0 --fresh-reference)
B2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --time=03:00:00 --output=$O/logs/boot_%A_%a.out --job-name=boot_v4 scripts/scmd_bootstrap.sbatch $O --cell $CELL --n 10000 --ntasks 272 --task-offset 240 --fresh-reference)
echo "bootstrap jobs $B1 $B2"
until [ "$(ls $O/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 272 ]; do sleep 30; done
sleep 30; scancel $B1 $B2 2>/dev/null
echo "bootstrap done $(date '+%F %T')"
cat $O/bootstrap/cell.json; echo
$PY -m isochronetoolbox.scmd bootstrap-collect --out $O 2>&1 | grep -v Warn
$PY $B/tuc_statistics.py $O 2>&1 | grep -v Warn
$PY $B/best3_and_null.py $O v4 2>&1 | grep -v Warn
echo "all done $(date '+%F %T')"
