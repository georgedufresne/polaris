#!/bin/bash
# After the 500-run grid: status, collect, the original posterior analysis on new results,
# comparison with OG (same 500 runs), plots of the 10 best (distinct runs).
B=/dartfs-hpc/scratch/f007058/scmd_tb; O=$B/og500
source $(conda info --base)/etc/profile.d/conda.sh; conda activate fidanka-env
cd /dartfs-hpc/rc/home/8/f007058/isochrone-toolbox
python -m isochronetoolbox.scmd status  --out $O --runs $O/runs.txt
python -m isochronetoolbox.scmd collect --out $O
A=$B/analysis_new500; mkdir -p $A/Bootstrap $A/out
cp $O/filtered_ks_scmd.csv $A/; cp $B/analysis_og500/initial_vars.csv $A/; cp $B/analysis_og500/Bootstrap/bootstrap_2d_ks_results.csv $A/Bootstrap/
SCMD_HOME=$A SCMD_OUTPUT_DIR=$A/out python /dartfs-hpc/scratch/f007058/scmd_noco_pipeline/analyze_scmd_only.py > $A/out/log.txt 2>&1 && tail -10 $A/out/log.txt
python $B/compare_new_vs_og.py
python -m isochronetoolbox.scmd plot --out $O --top 10 --distinct-runs
