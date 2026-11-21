#!/bin/bash
# tuc_final: the 47 Tuc DSEP4 grid with the settings from tuc_errchain: zero-point map (zp_ref 17.9-19.0,
# F606W share 0.08 -- nearly an F814W-only pattern), radii from the photometry before the ridge cut,
# PDMF -0.45 and binaries 0.009 (literature/scmd_cluster_inputs.csv), mu 13.15-13.55, keep 12.
# Waits for m15_final to finish (CPU limit), then grid (physics 240 + free 32), collect, bootstrap
# (10,000, fresh references, resampled ASTs), statistics, best-3 figure.
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export SCMD_PYTHON=$PY HOME=/dartfs-hpc/scratch/f007058 OMP_NUM_THREADS=1
B=/dartfs-hpc/scratch/f007058/scmd_tb
O=$B/tuc_final
export TUC_AGE_LIMITS=$O/tuc_age_limits_d4.json
cd $O/code
st() { $PY -m isochronetoolbox.scmd status --out $O --runs $O/runs.txt 2>/dev/null | tail -1; }
until grep -q "all done" $B/m15_final/finish.log 2>/dev/null || ! pgrep -f finish_m15_final.sh >/dev/null; do sleep 300; done
until ! squeue -h -u f007058 -n m15f_scmd,boot_m15f,tuc_chain,tuc_chain6,tuc_null | grep -q .; do sleep 120; done
echo "M15 final finished; starting 47 Tuc $(date '+%F %T')"
for r in $(cat $O/runs.txt); do [ -s $B/tuc_iso_d4v2/$r/isos/isochrones_example.txt ] || echo "missing isochrones: $r"; done
G1=$(sbatch --parsable --account=physics --array=0-239 --job-name=tucf_scmd --output=$O/logs/%A_%a.out scripts/scmd_worker.sbatch $O $O/runs.txt --max-hours 23)
G2=$(sbatch --parsable --account=free --array=0-31 --job-name=tucf_scmd --output=$O/logs/%A_%a.out scripts/scmd_worker.sbatch $O $O/runs.txt --max-hours 23)
echo "grid jobs $G1 $G2 $(date '+%F %T')" | tee -a $O/jobids.txt
until st | grep -q "'claimed': 0, 'todo': 0" || ! squeue -h -u f007058 -n tucf_scmd | grep -q .; do sleep 120; done
sleep 60; scancel -n tucf_scmd -u f007058 2>/dev/null; sleep 20
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
B1=$(sbatch --parsable --export=ALL --account=physics --array=0-239 --time=04:00:00 --output=$O/logs/boot_%A_%a.out --job-name=boot_tucf scripts/scmd_bootstrap.sbatch $O --cell $CELL --n 10000 --ntasks 272 --task-offset 0 --fresh-reference)
B2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --time=04:00:00 --output=$O/logs/boot_%A_%a.out --job-name=boot_tucf scripts/scmd_bootstrap.sbatch $O --cell $CELL --n 10000 --ntasks 272 --task-offset 240 --fresh-reference)
echo "bootstrap jobs $B1 $B2" | tee -a $O/jobids.txt
until [ "$(ls $O/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 272 ]; do sleep 30; done
sleep 30; scancel $B1 $B2 2>/dev/null
echo "bootstrap done $(date '+%F %T')"
cat $O/bootstrap/cell.json; echo
$PY -m isochronetoolbox.scmd bootstrap-collect --out $O 2>&1 | grep -v Warn
$PY $B/tucf_statistics.py $O 2>&1 | grep -v Warn
$PY $B/best3_and_null.py $O DSEP4 2>&1 | grep -v Warn
echo "all done $(date '+%F %T')"
