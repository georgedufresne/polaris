#!/bin/bash
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export SCMD_PYTHON=$PY HOME=/dartfs-hpc/scratch/f007058
O=/dartfs-hpc/scratch/f007058/scmd_tb/hugs/run10005
cd $O/code
until [ "$($PY -m isochronetoolbox.scmd status --out $O --runs $O/runs.txt 2>/dev/null | grep -o "'done': [0-9]*" | grep -o "[0-9]*$")" = "61" ]; do sleep 30; done
sleep 20; scancel -n hugs10005 -u f007058 2>/dev/null; sleep 10
$PY -m isochronetoolbox.scmd collect --out $O 2>&1 | grep -v Warn
J=$(sbatch --parsable --export=ALL --account=physics --array=0-60 --time=02:00:00 --output=$O/logs/boot_%A_%a.out --job-name=hugsboot scripts/scmd_bootstrap.sbatch $O --n 2000 --ntasks 61 --task-offset 0)
echo "bootstrap $J"
until [ "$(ls $O/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 61 ]; do sleep 20; done
sleep 20; scancel $J 2>/dev/null
cat $O/bootstrap/cell.json; echo
$PY -m isochronetoolbox.scmd bootstrap-collect --out $O 2>&1 | grep -v Warn
