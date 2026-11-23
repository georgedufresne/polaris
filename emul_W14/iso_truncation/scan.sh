#!/bin/bash
# Count isochrone points per age in every run of a W isochrone set; flag ages with < 90 % of the mean count of
# their two neighbours (truncated curves, e.g. Isochrones_W14/20374 at 7.2 Gyr: 926 vs ~2320 points).
# usage (array): scan.sh <iso_dir> <label>; task i takes runs i, i+n, ... (n = SLURM_ARRAY_TASK_COUNT)
#SBATCH --time=02:00:00 --cpus-per-task=1 --mem=2G
D=${1:?}; L=${2:?}; O=/dartfs-hpc/scratch/f007058/emul_W14/iso_truncation/output
n=${SLURM_ARRAY_TASK_COUNT:?}; i=${SLURM_ARRAY_TASK_ID:?}
ls $D | sort -n | awk -v n=$n -v i=$i '(NR-1) % n == i' | while read r; do
  f=$D/$r/isos/isochrones_example.txt; [ -s $f ] || { echo -e "$r\tmissing"; continue; }
  awk '!/^#/ && NF > 5 {c[$2 + 0]++} END {
       nk = 0; for (a in c) ages[++nk] = a + 0; asort(ages)
       bad = ""; for (j = 1; j <= nk; j++) {
         lo = (j > 1) ? c[ages[j-1]] : c[ages[j+1]]; hi = (j < nk) ? c[ages[j+1]] : c[ages[j-1]]; ref = (lo + hi) / 2
         if (c[ages[j]] < 0.9 * ref) bad = bad sprintf(" %.1f:%d/%d", ages[j] / 1e9, c[ages[j]], ref) }
       printf "%d ages\t%s\n", nk, (bad == "" ? "ok" : bad)}' $f | sed "s/^/$r\t/"
done > $O/${L}_$i.tsv
echo "done $(date '+%F %T')"
