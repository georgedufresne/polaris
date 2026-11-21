#!/bin/bash
# tuc_grid_d4: the v4 model (D: Treasury zero points, no short/long split, paired nearest-AST
# errors, setup.json identical to tuc_grid_v4's but for the isochrones) on the 1000 DSEP4 47 Tuc
# runs, isochrones tuc_iso_d4v2 (pyeep stopped at the RGB tip, tuc_d4v2_build/).
#  1. when the isochrone build (job BUILD) is done: grid workers on the runs that have isochrones
#     (physics 197 + free 32, so with the 43 DSEP4 reruns on physics the account stays at 240);
#  2. when the DSEP4 rerun (job RERUN, which builds its own isochrones) is done: 43 more workers
#     on those runs;
#  3. when every unit is done: collect, age limits, bootstrap (fresh references, resampled ASTs)
#     at the best cell within isochrone coverage, statistics, best-3 figure -- as finish_v4.sh.
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export SCMD_PYTHON=$PY HOME=/dartfs-hpc/scratch/f007058 OMP_NUM_THREADS=1
B=/dartfs-hpc/scratch/f007058/scmd_tb
O=$B/tuc_grid_d4
IB=$B/tuc_d4v2_build
BUILD=10731705
RERUN=10731704
cd $O/code
running() { squeue -h -j "$1" 2>/dev/null | grep -q .; }
have_iso() { for r in $(cat "$1"); do [ -s $B/tuc_iso_d4v2/$r/isos/isochrones_example.txt ] && echo $r; done; }
st() { $PY -m isochronetoolbox.scmd status --out $O --runs $O/runs.txt 2>/dev/null | tail -1; }

while running $BUILD; do sleep 60; done
echo "isochrone build done $(date '+%F %T')"
have_iso $IB/runs_complete.txt > $O/runs.txt
echo "runs with isochrones: $(wc -l < $O/runs.txt) of $(wc -l < $IB/runs_complete.txt)"
G1=$(sbatch --parsable --account=physics --array=0-196 --job-name=tuc_d4 --output=$O/logs/%A_%a.out scripts/scmd_worker.sbatch $O $O/runs.txt --max-hours 23)
G2=$(sbatch --parsable --account=free --array=0-31 --job-name=tuc_d4 --output=$O/logs/%A_%a.out scripts/scmd_worker.sbatch $O $O/runs.txt --max-hours 23)
echo "grid jobs $G1 $G2" | tee -a $O/jobids.txt

while running $RERUN; do sleep 120; done
echo "DSEP4 rerun done $(date '+%F %T')"
have_iso $IB/runs_rerun.txt > $O/runs_rerun.txt
n_new=$(wc -l < $O/runs_rerun.txt)
echo "rerun runs with isochrones: $n_new of $(wc -l < $IB/runs_rerun.txt)"
if [ "$n_new" -gt 0 ]; then
  cat $O/runs_rerun.txt >> $O/runs.txt
  G3=$(sbatch --parsable --account=physics --array=0-$((n_new - 1)) --job-name=tuc_d4 --output=$O/logs/%A_%a.out scripts/scmd_worker.sbatch $O $O/runs_rerun.txt --max-hours 23)
  echo "grid job (rerun runs) $G3" | tee -a $O/jobids.txt
fi

# done when nothing is claimed or left; stop waiting if no worker is left (a failed unit keeps its claim)
until st | grep -q "'claimed': 0, 'todo': 0" || ! squeue -h -u f007058 -n tuc_d4 | grep -q .; do sleep 120; done
sleep 60; scancel -n tuc_d4 -u f007058 2>/dev/null; sleep 20
echo "grid done $(date '+%F %T')"; st
$PY -m isochronetoolbox.scmd collect --out $O 2>&1 | grep -v Warn
$PY $IB/age_limits_d4.py $O/runs.txt $O/tuc_age_limits_d4.json
export TUC_AGE_LIMITS=$O/tuc_age_limits_d4.json
CELL=$($PY - <<PYEOF
import json, pandas as pd
lim = {int(k): v for k, v in json.load(open("$TUC_AGE_LIMITS")).items()}
g = pd.read_csv("$O/combined_2d_ks_results.csv")
g = g[g.age_gyr < g.simulation_id.map(lim) - 1e-6]
r = g.loc[g.ks_statistic.idxmin()]
print(f"{int(r.simulation_id)}:{r.age_gyr:.4f}:{r.distance_modulus:.4f}:{r.extinction_ebv:.4f}")
PYEOF
)
echo "bootstrap cell $CELL"
B1=$(sbatch --parsable --export=ALL --account=physics --array=0-239 --time=03:00:00 --output=$O/logs/boot_%A_%a.out --job-name=boot_d4 scripts/scmd_bootstrap.sbatch $O --cell $CELL --n 10000 --ntasks 272 --task-offset 0 --fresh-reference)
B2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --time=03:00:00 --output=$O/logs/boot_%A_%a.out --job-name=boot_d4 scripts/scmd_bootstrap.sbatch $O --cell $CELL --n 10000 --ntasks 272 --task-offset 240 --fresh-reference)
echo "bootstrap jobs $B1 $B2"
until [ "$(ls $O/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 272 ]; do sleep 30; done
sleep 30; scancel $B1 $B2 2>/dev/null
echo "bootstrap done $(date '+%F %T')"
cat $O/bootstrap/cell.json; echo
$PY -m isochronetoolbox.scmd bootstrap-collect --out $O 2>&1 | grep -v Warn
$PY $B/tuc_statistics.py $O 2>&1 | grep -v Warn
$PY $B/best3_and_null.py $O DSEP4 2>&1 | grep -v Warn
echo "all done $(date '+%F %T')"
