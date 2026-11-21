#!/bin/bash
# W 500-run grid: collect, posterior with the original analysis, W vs old-pipeline W (same runs, should agree),
# W vs OG (same runs: the emulator's shift), plots of the 10 best W fits (scored catalogues).
B=/dartfs-hpc/scratch/f007058/scmd_tb; W=$B/w500; R=/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox
source $(conda info --base)/etc/profile.d/conda.sh; conda activate fidanka-env
cd $R
python -m isochronetoolbox.scmd status  --out $W --runs $W/runs.txt
python -m isochronetoolbox.scmd collect --out $W
A=$B/analysis_w500; mkdir -p $A/Bootstrap $A/out
cp $W/filtered_ks_scmd.csv $A/; cp $B/analysis_og500/initial_vars.csv $A/; cp $B/analysis_og500/Bootstrap/bootstrap_2d_ks_results.csv $A/Bootstrap/
SCMD_HOME=$A SCMD_OUTPUT_DIR=$A/out python /dartfs-hpc/scratch/f007058/scmd_noco_pipeline/analyze_scmd_only.py > $A/out/log.txt 2>&1 && tail -9 $A/out/log.txt
echo; echo "=== W (toolbox) vs W (old scripts), same runs: should agree ==="
python tests/scmd/validate_against_reference.py --new $W --ref /dartfs-hpc/scratch/f007058/2d_ks_W --runs $W/runs.txt --report $B/w500_vs_oldW.txt
echo; echo "=== W (toolbox) vs OG DSEP isochrones (old tables), same runs: the emulator's shift ==="
python tests/scmd/validate_against_reference.py --new $W --ref /dartfs-hpc/rc/lab/C/ChaboyerB/2d_ks --runs $W/runs.txt --max-z 0 --report $B/w500_vs_og.txt
echo; echo "=== W vs OG, both run with the toolbox ==="
python tests/scmd/validate_against_reference.py --new $W --ref $B/og500 --runs $W/runs.txt --max-z 0 --report $B/w500_vs_og500toolbox.txt
python -m isochronetoolbox.scmd plot --out $W --top 10 --distinct-runs
