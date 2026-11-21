#!/bin/bash
# Copy the notes written on scratch into ~/notes/scratch (same layout). Home is not purged and has weekly snapshots.
# Notes = *.md, README*, NOTE*, and text summaries (*.txt <= 200 KB; machine lists such as run or job-id lists excluded).
# Never deletes anything in ~/notes. Stops if home has less than 1 GB free.
set -euo pipefail
SRC=/dartfs-hpc/scratch/f007058/
DST=/dartfs-hpc/rc/home/8/f007058/notes/scratch/
free_gb=$(df -BG --output=avail /dartfs-hpc/rc/home/8/f007058 | tail -1 | tr -dc 0-9)
[ "$free_gb" -ge 1 ] || { echo "home has ${free_gb} GB free; notes not synced"; exit 1; }
mkdir -p "$DST"
rsync -rlt -m --max-size=200K \
  --exclude='.git/' --exclude='__pycache__/' \
  --exclude='/envs/' --exclude='/.cache/' --exclude='/.config/' --exclude='/isochrone-toolbox/' \
  --exclude='/repository-instructions-main/' --exclude='/Files_from_Martin/' \
  --exclude='/sCMD_*/' --exclude='/Isochrones_*/' --exclude='/Tracks_*/' --exclude='/plots_for_2d_*/' --exclude='/47Tuc_iso/' \
  --exclude='/2d_ks_*/' --exclude='/emul_W14/CONF1DENCE-W_*/' --exclude='/M15/zenodo_ying2025/' --exclude='/M15/iso_*/' \
  --exclude='/M15/tracks_*/' --exclude='/tuc_dsep4_batch/work/' --exclude='/M15/dsep4/work/' --exclude='/M15/dsep4/env/DSEPmisc/' \
  --exclude='/scmd_tb/tuc_iso*/[0-9]*/' --exclude='/scmd_tb/tuc_iso*/w[0-9]*/' --exclude='/audit_ngc6791/pyeep_work/' \
  --exclude='units/' --exclude='claims/' --exclude='logs/' --exclude='logs_*/' --exclude='pyeep/' --exclude='code/' --exclude='tests/' \
  --exclude='jobids*' --exclude='*runs*.txt' --exclude='jobmap.txt' --exclude='*cancelled*' --exclude='input*.txt' \
  --exclude='isochrones*.txt' --exclude='requirements*.txt' \
  --include='*/' --include='*.md' --include='README*' --include='NOTE*' --include='*.txt' --exclude='*' \
  "$SRC" "$DST"
echo "notes synced $(date '+%F %T'): $(find "$DST" -type f | wc -l) files, $(du -sh "$DST" | cut -f1) in $DST"
