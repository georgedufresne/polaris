#!/bin/bash
# Copy the code written on Polaris scratch (/dartfs-hpc/scratch/f007058) into this repo, keeping the scratch layout
# (so paths in the scripts map 1:1). Re-run any time, then commit. Never deletes anything here.
# Copied: *.py, shell/Slurm drivers, READMEs/notes (*.md), small configs (json/toml/nml/yaml), each <= 1 MB.
# Not copied: data, isochrones, tracks, sCMDs, plots, logs; other people's code (Files_from_Martin, the CONF1DENCE-W
# snapshot, Ying's Zenodo files, *verbatim* copies); separate repos (isochrone-toolbox); envs and caches.
set -euo pipefail
SRC=/dartfs-hpc/scratch/f007058/
DST=$(cd "$(dirname "$0")" && pwd)/
rsync -rlt -m --max-size=1M \
  --exclude='.git/' --exclude='__pycache__/' --exclude='.ipynb_checkpoints/' \
  --exclude='/envs/' --exclude='/.cache/' --exclude='/.config/' --exclude='/isochrone-toolbox/' \
  --exclude='/repository-instructions-main/' --exclude='/Files_from_Martin/' --exclude='/storage_survey/' \
  --exclude='/sCMD_*/' --exclude='/Isochrones_*/' --exclude='/Tracks_*/' --exclude='/plots_for_2d_*/' --exclude='/47Tuc_iso/' \
  --exclude='/emul_W14/CONF1DENCE-W_*/' --exclude='/M15/zenodo_ying2025/' --exclude='/M15/iso_*/' --exclude='/M15/tracks_*/' \
  --exclude='*verbatim*' --exclude='units/' --exclude='claims/' --exclude='logs/' --exclude='logs_*/' \
  --exclude='/tuc_dsep4_batch/work/' --exclude='/M15/dsep4/work/' --exclude='/M15/dsep4/env/DSEPmisc/' \
  --exclude='/scmd_tb/tuc_iso*/[0-9]*/' --exclude='/scmd_tb/tuc_iso*/w[0-9]*/' \
  --include='*/' \
  --include='*.py' --include='*.sh' --include='*.batch' --include='*.sbatch' --include='*.slurm' \
  --include='*.md' --include='*.json' --include='*.toml' --include='*.nml' --include='*.yaml' --include='*.yml' \
  --exclude='*' \
  "$SRC" "$DST"
echo "synced $(date '+%F %T'): $(find "$DST" -type f -not -path '*/.git/*' | wc -l) files, $(du -sh --exclude=.git "$DST" | cut -f1)"
