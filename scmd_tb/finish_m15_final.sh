#!/bin/bash
# m15_final: the M15 DSEP4 grid with the settings from m15_errchain (chain_summary.csv):
#   zero-point map on both bands (zp_ref 19.3-20.6, F606W share 0.57, diag/diag3.py), radii from the
#   photometry before the ridge cut, PDMF -1.00 and binaries 0.017 (literature/scmd_cluster_inputs.csv),
#   E 0.050-0.120 (the E(V-I) 0.08-0.15 prior converted to F606W-F814W, and below), keep 12.
#  1. grid workers: physics 200 + free 32 (40 physics CPUs left for the 47 Tuc tests; more added later)
#  2. collect, bootstrap (10,000, fresh references, resampled ASTs) at the lowest-KS cell, statistics
#     (all E, and E >= 0.08 only), best-3 figure.
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export SCMD_PYTHON=$PY HOME=/dartfs-hpc/scratch/f007058 OMP_NUM_THREADS=1
B=/dartfs-hpc/scratch/f007058/scmd_tb
O=$B/m15_final
export TUC_AGE_LIMITS=/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json
cd $O/code
st() { $PY -m isochronetoolbox.scmd status --out $O --runs $O/runs.txt 2>/dev/null | tail -1; }
G1=$(sbatch --parsable --account=physics --array=0-199 --job-name=m15f_scmd --output=$O/logs/%A_%a.out scripts/scmd_worker.sbatch $O $O/runs.txt --max-hours 23)
G2=$(sbatch --parsable --account=free --array=0-31 --job-name=m15f_scmd --output=$O/logs/%A_%a.out scripts/scmd_worker.sbatch $O $O/runs.txt --max-hours 23)
echo "grid jobs $G1 $G2 $(date '+%F %T')" | tee -a $O/jobids.txt
until st | grep -q "'claimed': 0, 'todo': 0" || ! squeue -h -u f007058 -n m15f_scmd | grep -q .; do sleep 120; done
sleep 60; scancel -n m15f_scmd -u f007058 2>/dev/null; sleep 20
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
B1=$(sbatch --parsable --export=ALL --account=physics --array=0-239 --time=03:00:00 --output=$O/logs/boot_%A_%a.out --job-name=boot_m15f scripts/scmd_bootstrap.sbatch $O --cell $CELL --n 10000 --ntasks 272 --task-offset 0 --fresh-reference)
B2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --time=03:00:00 --output=$O/logs/boot_%A_%a.out --job-name=boot_m15f scripts/scmd_bootstrap.sbatch $O --cell $CELL --n 10000 --ntasks 272 --task-offset 240 --fresh-reference)
echo "bootstrap jobs $B1 $B2" | tee -a $O/jobids.txt
until [ "$(ls $O/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 272 ]; do sleep 30; done
sleep 30; scancel $B1 $B2 2>/dev/null
echo "bootstrap done $(date '+%F %T')"
cat $O/bootstrap/cell.json; echo
$PY -m isochronetoolbox.scmd bootstrap-collect --out $O 2>&1 | grep -v Warn
$PY $B/m15f_statistics.py $O 2>&1 | grep -v Warn
echo "=== E >= 0.08 only"
EVI_MIN=0.08 POST_SUFFIX=_Efloor $PY $B/m15f_statistics.py $O 2>&1 | grep -v Warn
$PY $B/m15_best3_and_null.py $O DSEP4 2>&1 | grep -v Warn
echo "all done $(date '+%F %T')"
