#!/bin/bash
source $(conda info --base)/etc/profile.d/conda.sh; conda activate fidanka-env
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
H=/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/orientation
cd $H
test -f $H/ref_split.npy || { echo "no reference"; exit 1; }
python run_part.py null --ids 0:40:4 --out $H/null_0.jsonl --save 1 --validate 1 > $H/null_0.err 2>&1 &
python run_part.py null --ids 1:40:4 --out $H/null_1.jsonl --save 1 --validate 0 > $H/null_1.err 2>&1 &
python run_part.py null --ids 2:40:4 --out $H/null_2.jsonl --save 0 --validate 0 > $H/null_2.err 2>&1 &
python run_part.py null --ids 3:40:4 --out $H/null_3.jsonl --save 0 --validate 0 > $H/null_3.err 2>&1 &
wait
echo PHASE2 DONE
