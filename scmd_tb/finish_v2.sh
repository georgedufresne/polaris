#!/bin/bash
# tuc_grid_v2: collect, bootstrap (fresh references, resampled ASTs) at the best cell, statistics.
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export SCMD_PYTHON=$PY HOME=/dartfs-hpc/scratch/f007058
O=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_v2
cd $O/code
until [ "$($PY -m isochronetoolbox.scmd status --out $O --runs $O/runs.txt 2>/dev/null | grep -o "'done': [0-9]*" | grep -o "[0-9]*$")" = "679" ]; do sleep 60; done
sleep 60; scancel -n tuc_v2 -u f007058 2>/dev/null; sleep 20
$PY -m isochronetoolbox.scmd collect --out $O 2>&1 | grep -v Warn
B1=$(sbatch --parsable --export=ALL --account=physics --array=0-239 --time=03:00:00 --output=$O/logs/boot_%A_%a.out --job-name=boot_v2 scripts/scmd_bootstrap.sbatch $O --n 10000 --ntasks 272 --task-offset 0 --fresh-reference)
B2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --time=03:00:00 --output=$O/logs/boot_%A_%a.out --job-name=boot_v2 scripts/scmd_bootstrap.sbatch $O --n 10000 --ntasks 272 --task-offset 240 --fresh-reference)
echo "bootstrap jobs $B1 $B2"
until [ "$(ls $O/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 272 ]; do sleep 30; done
sleep 30; scancel $B1 $B2 2>/dev/null
cat $O/bootstrap/cell.json; echo
$PY -m isochronetoolbox.scmd bootstrap-collect --out $O 2>&1 | grep -v Warn
$PY /dartfs-hpc/scratch/f007058/scmd_tb/tuc_statistics.py $O 2>&1 | grep -v Warn
