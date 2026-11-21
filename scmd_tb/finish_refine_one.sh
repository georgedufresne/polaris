#!/bin/bash
# Steps 3-5 of finish_refine.sh for one cluster (m15 | tuc), so M15 does not wait for the 47 Tuc refinement.
# finish_refine.sh was stopped on 2026-11-05 in step 2, after the M15 refinement had finished.
#  tuc: waits for the refinement to finish first, then validates.
#  then: merge -> <grid>/combined_2d_ks_results_refined.csv; <grid>_refined/ (setup with the E override,
#  code link, refined rows); 10,000 bootstraps at the lowest-KS cell below the age limits (136 tasks:
#  physics 120 + free 16, so both clusters together stay within 240 + 32); statistics; best-3 figure.
X=${1:?usage: finish_refine_one.sh m15|tuc}
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export SCMD_PYTHON=$PY HOME=/dartfs-hpc/scratch/f007058 OMP_NUM_THREADS=1
B=/dartfs-hpc/scratch/f007058/scmd_tb
K=$B/ks_refine
G=$B/${X}_final
V=$B/${X}_final_refined
case $X in
    m15) LIM=/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json ;;
    tuc) LIM=$G/tuc_age_limits_d4.json ;;
    *) echo "unknown cluster $X"; exit 1 ;;
esac
echo $$ > $K/driver_$X.pid

rs() { $PY $K/refine.py $G status 2>/dev/null | tail -1; }
if [ $X = tuc ]; then
    echo "=== $X: waiting for the refinement $(date '+%F %T')"
    until rs | grep -q "'claimed': 0, 'todo': 0" || ! squeue -h -u f007058 -n ks_refine | grep -q .; do sleep 120; done
    sleep 30; scancel -n ks_refine -u f007058 2>/dev/null
    echo "refinement done $(date '+%F %T')"; rs
fi

echo "=== $X: validation and merge $(date '+%F %T')"
$PY $K/validate.py $X
$PY $K/refine.py $G merge || exit 1

echo "=== $X: bootstrap $(date '+%F %T')"
mkdir -p $V/logs $V/analysis
cp $G/setup.json $G/runs.txt $V/
if [ -f $G/refine/setup_override.json ]; then
    $PY -c "import json; s = json.load(open('$V/setup.json')); s.update(json.load(open('$G/refine/setup_override.json'))); json.dump(s, open('$V/setup.json', 'w'), indent=1)"
    echo "setup override: $(cat $G/refine/setup_override.json)"
fi
ln -sfn $G/code $V/code
cp $G/combined_2d_ks_results_refined.csv $V/combined_2d_ks_results.csv
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
cd $G/code
B1=$(sbatch --parsable --export=ALL --account=physics --array=0-119 --time=02:00:00 --output=$V/logs/boot_%A_%a.out --job-name=boot_$X scripts/scmd_bootstrap.sbatch $V --cell $CELL --n 10000 --ntasks 136 --task-offset 0 --fresh-reference)
B2=$(sbatch --parsable --export=ALL --account=free --array=0-15 --time=02:00:00 --output=$V/logs/boot_%A_%a.out --job-name=boot_$X scripts/scmd_bootstrap.sbatch $V --cell $CELL --n 10000 --ntasks 136 --task-offset 120 --fresh-reference)
echo "bootstrap jobs $B1 $B2" | tee -a $V/jobids.txt $G/jobids.txt
until [ "$(ls $V/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 136 ]; do
    squeue -h -u f007058 -n boot_$X | grep -q . || break
    sleep 60
done
sleep 30; scancel -n boot_$X -u f007058 2>/dev/null
echo "bootstrap done $(date '+%F %T'): $(ls $V/bootstrap/part_*.csv 2>/dev/null | wc -l) parts"

echo "=== $X: statistics $(date '+%F %T')"
$PY -m isochronetoolbox.scmd bootstrap-collect --out $V 2>&1 | grep -v Warn
if [ $X = m15 ]; then
    echo "--- M15, all E (0.000-0.120)"
    $PY $B/m15f_statistics.py $V 2>&1 | grep -v Warn
    echo "--- M15, E >= 0.05 (the original grid floor)"
    EVI_MIN=0.05 POST_SUFFIX=_E05 $PY $B/m15f_statistics.py $V 2>&1 | grep -v Warn
    echo "--- M15, E >= 0.08 (literature E(B-V) 0.10 +/- 0.02)"
    EVI_MIN=0.08 POST_SUFFIX=_Elit $PY $B/m15f_statistics.py $V 2>&1 | grep -v Warn
    TUC_AGE_LIMITS=$LIM $PY $B/m15_best3_and_null.py $V DSEP4 2>&1 | grep -v Warn
else
    $PY $B/tucf_statistics.py $V 2>&1 | grep -v Warn
    TUC_AGE_LIMITS=$LIM $PY $B/best3_and_null.py $V DSEP4 2>&1 | grep -v Warn
fi
echo "all done $(date '+%F %T')"
