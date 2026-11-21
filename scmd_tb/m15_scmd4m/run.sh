#!/bin/bash
# m15_scmd4m: M15 rescored with 4,000,000-star model sCMDs (as Ying et al. 2025), for two baselines:
#   prod/  production (m15_final_refined): zero-point map on, E 0.000-0.120
#   nozp/  m15_yingbox_nozp_refined: Ying's E box 0.080-0.150, no zero-point map
# The only change from each baseline: n_model 4,000,000 (ridge still traced on an observed-sized sample) and the
# model ECDF counted exactly (ks_model_ecdf "count"). Waits for ../tuc_yingbox_nozp/run.sh (CPU limit), then:
#  1. refine4m.py workers (physics 240 + free 32): 8-neighbour KS descent with 4M sCMDs from every 32k
#     (run, age) minimum with KS <= 2.5 x the grid best; nozp first, then prod; merge
#  2. per baseline: 10,000-realisation bootstrap at the lowest 4M cell (fresh N-star references with all ASTs,
#     4M-star models with resampled ASTs); statistics; comparison figure
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export SCMD_PYTHON=$PY HOME=/dartfs-hpc/scratch/f007058 OMP_NUM_THREADS=1
B=/dartfs-hpc/scratch/f007058/scmd_tb
H=$B/m15_scmd4m
T=$B/tuc_yingbox_nozp
LIM=/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json
export TUC_AGE_LIMITS=$LIM
echo $$ > $H/driver.pid

echo "=== waiting for tuc_yingbox_nozp $(date '+%F %T')"
until grep -q "all done" $T/run.log 2>/dev/null || ! kill -0 "$(cat $T/driver.pid 2>/dev/null)" 2>/dev/null; do sleep 300; done
until ! squeue -h -u f007058 -n tucyb_scmd,tucyb_refine,boot_tucyb | grep -q .; do sleep 120; done
echo "47 Tuc finished; starting $(date '+%F %T')"

echo "=== 4M descent $(date '+%F %T')"
st() { for X in nozp prod; do echo "$X $($PY $H/refine4m.py $H/$X status 2>/dev/null | tail -1)"; done; }
W1=$(sbatch --parsable --account=physics --array=0-239 --output=$H/logs/w_%A_%a.out $H/refine4m.sbatch $H/nozp $H/prod)
W2=$(sbatch --parsable --account=free --array=0-31 --output=$H/logs/w_%A_%a.out $H/refine4m.sbatch $H/nozp $H/prod)
echo "worker jobs $W1 $W2" | tee -a $H/jobids.txt
until [ "$(st | grep -c "'claimed': 0, 'todo': 0")" -eq 2 ] || ! squeue -h -u f007058 -n m15_4m | grep -q .; do sleep 120; done
sleep 30; scancel -n m15_4m -u f007058 2>/dev/null
echo "descent done $(date '+%F %T')"; st
for X in nozp prod; do echo "--- merge $X"; $PY $H/refine4m.py $H/$X merge || exit 1; done

cd $H/code
for X in nozp prod; do
    V=$H/$X
    echo "=== bootstrap $X $(date '+%F %T')"
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
    B1=$(sbatch --parsable --export=ALL --account=physics --array=0-239 --time=02:00:00 --output=$V/logs/boot_%A_%a.out --job-name=boot_m15_4m scripts/scmd_bootstrap.sbatch $V --cell $CELL --n 10000 --ntasks 272 --task-offset 0 --fresh-reference)
    B2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --time=02:00:00 --output=$V/logs/boot_%A_%a.out --job-name=boot_m15_4m scripts/scmd_bootstrap.sbatch $V --cell $CELL --n 10000 --ntasks 272 --task-offset 240 --fresh-reference)
    echo "bootstrap jobs $B1 $B2" | tee -a $H/jobids.txt
    until [ "$(ls $V/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 272 ]; do
        squeue -h -u f007058 -n boot_m15_4m | grep -q . || break
        sleep 60
    done
    sleep 30; scancel -n boot_m15_4m -u f007058 2>/dev/null
    echo "bootstrap done $(date '+%F %T'): $(ls $V/bootstrap/part_*.csv 2>/dev/null | wc -l) parts"
    $PY -m isochronetoolbox.scmd bootstrap-collect --out $V 2>&1 | grep -v Warn
done

echo "=== statistics $(date '+%F %T')"
echo "--- nozp: Ying's box, mu 15.30-15.60, E 0.080-0.150 (the whole grid)"
$PY $B/m15f_statistics.py $H/nozp 2>&1 | grep -v Warn
echo "--- prod: all E (0.000-0.120)"
$PY $B/m15f_statistics.py $H/prod 2>&1 | grep -v Warn
echo "--- prod: box mu 15.3-15.6, E 0.08-0.15"
MU_MIN=15.3 MU_MAX=15.6 EVI_MIN=0.08 EVI_MAX=0.15 POST_SUFFIX=_DB $PY $B/m15f_statistics.py $H/prod 2>&1 | grep -v Warn
echo "--- prod: box mu 15.3-15.6, E 0.062-0.117 (x 0.78)"
MU_MIN=15.3 MU_MAX=15.6 EVI_MIN=0.062 EVI_MAX=0.117 POST_SUFFIX=_DBconv $PY $B/m15f_statistics.py $H/prod 2>&1 | grep -v Warn
$PY $H/compare_ks_vs_age.py 2>&1 | grep -v Warn
for X in nozp prod; do $PY $B/refined_plots.py $H/$X $LIM "M15 4M sCMDs ($X)" 2>&1 | grep -v Warn; done
echo "all done $(date '+%F %T')"
