#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh; conda activate fidanka-env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
H=/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/orientation
cd $H
python run_part.py ref  --out $H/ref.jsonl                 > $H/ref.err 2>&1 &
python run_part.py grid --out $H/grid.jsonl --save 1       > $H/grid.err 2>&1 &
python run_part.py data --ids 1000:1006 --out $H/data_a.jsonl --save 1 > $H/data_a.err 2>&1 &
python run_part.py data --ids 1006:1012 --out $H/data_b.jsonl --save 1 --validate 0 > $H/data_b.err 2>&1 &
wait
echo PHASE1 DONE
