#!/bin/bash
# Fresh 500-run sCMD grid with the Kerberos-free setup (everything on scratch).
set -e
O=/dartfs-hpc/scratch/f007058/scmd_tb/og500
R=/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox
rm -rf $O/units $O/claims $O/2d_ks $O/logs $O/plots $O/combined_2d_ks_results.csv $O/filtered_ks_scmd.csv $O/kde_js_scores.csv.gz
mkdir -p $O/logs
rsync -a --delete --exclude __pycache__ $R/isochronetoolbox $R/data $R/scripts $O/code/
export SCMD_PYTHON=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export HOME=/dartfs-hpc/scratch/f007058
cd $O/code
J1=$(sbatch --parsable --export=ALL --account=physics --array=0-239 --time=04:00:00 --output=$O/logs/%A_%a.out --job-name=scmd_og500 scripts/scmd_worker.sbatch $O $O/runs.txt)
J2=$(sbatch --parsable --export=ALL --account=free --array=0-31 --time=04:00:00 --output=$O/logs/%A_%a.out --job-name=scmd_og500 scripts/scmd_worker.sbatch $O $O/runs.txt)
echo "physics $J1  free $J2"
