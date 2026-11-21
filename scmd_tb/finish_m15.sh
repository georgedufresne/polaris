#!/bin/bash
# m15_grid_d4: the 47 Tuc DSEP4 model (tuc_grid_d4/setup.json: AST errors from the nearest
# artificial stars, Treasury zero points, binaries 2 %, alpha = -0.45, our ridge cut at 0.15 in
# normalised perpendicular distance, 95th-percentile 2D KS) on the 500 M15 DSEP4 runs
# (M15/iso_dsep4), with M15's inputs: M. Ying's cleaned M15.dat and artificial stars, SGB F606W
# 18.64, mu 15.30-15.60 (Bhardwaj et al. 2021, Feuillet et al. 2014), E(V-I) 0.08-0.15 (Dutra &
# Bica 2000), ages 8.0-15.0 Gyr (every run's isochrones are complete to 15.0).
#  1. grid workers: physics 240 + free 32;
#  2. when every unit is done: collect, bootstrap (10,000 fresh references, resampled ASTs) at
#     the lowest-KS cell, statistics, best-3 figure -- as finish_d4.sh.
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export SCMD_PYTHON=$PY HOME=/dartfs-hpc/scratch/f007058 OMP_NUM_THREADS=1
B=/dartfs-hpc/scratch/f007058/scmd_tb
O=$B/m15_grid_d4
export TUC_AGE_LIMITS=/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json
cd $O/code
st() { $PY -m isochronetoolbox.scmd status --out $O --runs $O/runs.txt 2>/dev/null | tail -1; }

for r in $(seq 10000 10499); do [ -s /dartfs-hpc/scratch/f007058/M15/iso_dsep4/$r/isos/isochrones_example.txt ] && echo $r; done > $O/runs.txt
echo "runs with isochrones: $(wc -l < $O/runs.txt)"
G1=$(sbatch --parsable --account=physics --array=0-239 --job-name=m15_scmd --output=$O/logs/%A_%a.out scripts/scmd_worker.sbatch $O $O/runs.txt --max-hours 23)
G2=$(sbatch --parsable --account=free --array=0-31 --job-name=m15_scmd --output=$O/logs/%A_%a.out scripts/scmd_worker.sbatch $O $O/runs.txt --max-hours 23)
echo "grid jobs $G1 $G2 $(date '+%F %T')" | tee -a $O/jobids.txt

# done when nothing is claimed or left; stop waiting if no worker is left (a failed unit keeps its claim)
until st | grep -q "'claimed': 0, 'todo': 0" || ! squeue -h -u f007058 -n m15_scmd | grep -q .; do sleep 120; done
sleep 60; scancel -n m15_scmd -u f007058 2>/dev/null; sleep 20
echo "grid done $(date '+%F %T')"; st
$PY -m isochronetoolbox.scmd collect --out $O 2>&1 | grep -v Warn
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
B1=$(sbatch --parsable --export=ALL --account=physics --array=0-239 --time=03:00:00 --output=$O/logs/boot_%A_%a.out --job-name=boot_m15 scripts/scmd_bootstrap.sbatch $O --cell $CELL --n 10000 --ntasks 272 --task-offset 0 --fresh-reference)
B2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --time=03:00:00 --output=$O/logs/boot_%A_%a.out --job-name=boot_m15 scripts/scmd_bootstrap.sbatch $O --cell $CELL --n 10000 --ntasks 272 --task-offset 240 --fresh-reference)
echo "bootstrap jobs $B1 $B2" | tee -a $O/jobids.txt
until [ "$(ls $O/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 272 ]; do sleep 30; done
sleep 30; scancel $B1 $B2 2>/dev/null
echo "bootstrap done $(date '+%F %T')"
cat $O/bootstrap/cell.json; echo
$PY -m isochronetoolbox.scmd bootstrap-collect --out $O 2>&1 | grep -v Warn
$PY $B/m15_statistics.py $O 2>&1 | grep -v Warn
$PY $B/m15_best3_and_null.py $O DSEP4 2>&1 | grep -v Warn
echo "all done $(date '+%F %T')"
