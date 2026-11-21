#!/bin/bash
# One-track check of the LMonteCarlo fix: run 10357 ([alpha/Fe] +0.2), 0.78 Msun, fixed namelists (env/nml_mc).
D=/dartfs-hpc/scratch/f007058/M15/dsep4
export HOME=/dartfs-hpc/scratch/f007058 LD_LIBRARY_PATH= NML=$D/env/nml_mc STOREBASE=$D/test_mcfix/tracks
cd $D
time ./run_m15_dsep4_one.sh 10357 -2.27 780 > $D/test_mcfix/10357.780.log 2>&1
echo "rc $?"
