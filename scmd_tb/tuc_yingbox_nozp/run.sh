#!/bin/bash
# tuc_yingbox_nozp: the production 47 Tuc pipeline (tuc_final -> tuc_final_refined) on the same 1000 runs with
# one change: no spatial zero-point map (zp_correct false). Production's E range, 0.000-0.050 applied directly
# as the F606W-F814W shift, is already M. Ying's box, so the M15 box change has no 47 Tuc counterpart; his mu
# range 13.15-13.45 is applied in the _D17 posteriors, as for production. Everything else is production.
# Waits for ../m15_yingbox_nozp/run.sh to finish (CPU limit), then:
#  1. grid: physics 240 + free 32 workers; collect
#  2. KS refinement (E 0.000-0.050 in 0.001); merge
#  3. tuc_yingbox_nozp_refined/: setup with the E override, code link, refined rows; bootstrap (272 tasks)
#  4. statistics (whole grid; _D17 = Ying's box mu 13.15-13.45, E 0-0.05; _D17conv = E 0-0.039), figures
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export SCMD_PYTHON=$PY HOME=/dartfs-hpc/scratch/f007058 OMP_NUM_THREADS=1
B=/dartfs-hpc/scratch/f007058/scmd_tb
K=$B/ks_refine
O=$B/tuc_yingbox_nozp
V=$B/tuc_yingbox_nozp_refined
M=$B/m15_yingbox_nozp
LIM=$B/tuc_final/tuc_age_limits_d4.json
export TUC_AGE_LIMITS=$LIM
echo $$ > $O/driver.pid
cd $O/code

echo "=== waiting for m15_yingbox_nozp $(date '+%F %T')"
until grep -q "all done" $M/run.log 2>/dev/null || ! kill -0 "$(cat $M/driver.pid 2>/dev/null)" 2>/dev/null; do sleep 300; done
until ! squeue -h -u f007058 -n m15yb_scmd,m15yb_refine,boot_m15yb | grep -q .; do sleep 120; done
echo "M15 finished; starting 47 Tuc $(date '+%F %T')"
for r in $(cat $O/runs.txt); do [ -s $B/tuc_iso_d4v2/$r/isos/isochrones_example.txt ] || echo "missing isochrones: $r"; done

echo "=== grid $(date '+%F %T')"
st() { $PY -m isochronetoolbox.scmd status --out $O --runs $O/runs.txt 2>/dev/null | tail -1; }
G1=$(sbatch --parsable --account=physics --array=0-239 --job-name=tucyb_scmd --output=$O/logs/%A_%a.out scripts/scmd_worker.sbatch $O $O/runs.txt --max-hours 23)
G2=$(sbatch --parsable --account=free --array=0-31 --job-name=tucyb_scmd --output=$O/logs/%A_%a.out scripts/scmd_worker.sbatch $O $O/runs.txt --max-hours 23)
echo "grid jobs $G1 $G2 $(date '+%F %T')" | tee -a $O/jobids.txt
until st | grep -q "'claimed': 0, 'todo': 0" || ! squeue -h -u f007058 -n tucyb_scmd | grep -q .; do sleep 120; done
sleep 60; scancel -n tucyb_scmd -u f007058 2>/dev/null; sleep 20
echo "grid done $(date '+%F %T')"; st
$PY -m isochronetoolbox.scmd collect --out $O 2>&1 | grep -v Warn

echo "=== refinement $(date '+%F %T')"
mkdir -p $O/refine/logs $O/refine/claims
rs() { $PY $K/refine.py $O status 2>/dev/null | tail -1; }
R1=$(sbatch --parsable --export=ALL --account=physics --array=0-239 --job-name=tucyb_refine --output=$O/refine/logs/%A_%a.out $K/refine.sbatch $O $LIM)
R2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --job-name=tucyb_refine --output=$O/refine/logs/%A_%a.out $K/refine.sbatch $O $LIM)
echo "refine jobs $R1 $R2" | tee -a $O/jobids.txt
until rs | grep -q "'claimed': 0, 'todo': 0" || ! squeue -h -u f007058 -n tucyb_refine | grep -q .; do sleep 120; done
sleep 30; scancel -n tucyb_refine -u f007058 2>/dev/null
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
B1=$(sbatch --parsable --export=ALL --account=physics --array=0-239 --time=04:00:00 --output=$V/logs/boot_%A_%a.out --job-name=boot_tucyb scripts/scmd_bootstrap.sbatch $V --cell $CELL --n 10000 --ntasks 272 --task-offset 0 --fresh-reference)
B2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --time=04:00:00 --output=$V/logs/boot_%A_%a.out --job-name=boot_tucyb scripts/scmd_bootstrap.sbatch $V --cell $CELL --n 10000 --ntasks 272 --task-offset 240 --fresh-reference)
echo "bootstrap jobs $B1 $B2" | tee -a $V/jobids.txt $O/jobids.txt
until [ "$(ls $V/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 272 ]; do
    squeue -h -u f007058 -n boot_tucyb | grep -q . || break
    sleep 60
done
sleep 30; scancel -n boot_tucyb -u f007058 2>/dev/null
echo "bootstrap done $(date '+%F %T'): $(ls $V/bootstrap/part_*.csv 2>/dev/null | wc -l) parts"

echo "=== statistics $(date '+%F %T')"
$PY -m isochronetoolbox.scmd bootstrap-collect --out $V 2>&1 | grep -v Warn
echo "--- whole grid: mu 13.15-13.55, E 0.000-0.050"
$PY $B/tucf_statistics.py $V 2>&1 | grep -v Warn
echo "--- Ying's box: mu 13.15-13.45, E 0-0.05"
MU_MIN=13.15 MU_MAX=13.45 EVI_MIN=0 EVI_MAX=0.05 POST_SUFFIX=_D17 $PY $B/tucf_statistics.py $V 2>&1 | grep -v Warn
echo "--- box converted: mu 13.15-13.45, E 0-0.039"
MU_MIN=13.15 MU_MAX=13.45 EVI_MIN=0 EVI_MAX=0.039 POST_SUFFIX=_D17conv $PY $B/tucf_statistics.py $V 2>&1 | grep -v Warn
$PY $B/best3_and_null.py $V DSEP4 2>&1 | grep -v Warn
$PY $B/refined_plots.py $V $LIM "47 Tuc, no zero-point map" 2>&1 | grep -v Warn
echo "all done $(date '+%F %T')"
