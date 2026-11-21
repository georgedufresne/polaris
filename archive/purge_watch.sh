#!/bin/bash
#SBATCH --job-name=purge_watch
#SBATCH --account=free
#SBATCH --partition=standard
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --time=06:00:00
#SBATCH --output=/dartfs-hpc/rc/lab/C/ChaboyerB/george/polaris/archive/logs/purge_watch_%j.log
# What on scratch will the 45-day purge reach, and when? Read-only (find/stat do not touch the files).
# Age of a file = days since max(mtime, atime). VAST does not update atime on read, so in practice this is mtime.
# Report: george/purge_watch/report_<date>.txt and latest.txt; per top-level folder, GB and files that are
# already past 45 days, cross 45 days within 7 / 14 / 28 days, or are safer; plus the folder's status in kept.tsv.
set -uo pipefail
S=/dartfs-hpc/scratch/f007058
G=/dartfs-hpc/rc/lab/C/ChaboyerB/george
A=$G/polaris/archive
OUT=$G/purge_watch; mkdir -p $OUT
P=${SLURM_CPUS_PER_TASK:-8}
NOW=$(date +%s); DAY=$(date +%F)
RAW=$OUT/.raw_$DAY.tsv

# work units: every depth-2 path plus depth-1 files, so the huge folders are walked in parallel
{ find $S -mindepth 1 -maxdepth 1 -type f -printf '%P\n'; find $S -mindepth 2 -maxdepth 2 -printf '%P\n'; } |
  xargs -d '\n' -P $P -I{} sh -c '
    find "'$S'/{}" -type f -printf "%s %T@ %A@\n" 2>/dev/null |
    awk -v u="{}" -v now='$NOW' '"'"'
      { a = ($2 > $3 ? $2 : $3); d = (now - a) / 86400
        k = (d >= 45) ? 0 : (d >= 38) ? 1 : (d >= 31) ? 2 : (d >= 17) ? 3 : 4
        b[k] += $1; n[k]++; if (a > newest) newest = a }
      END { if (n[0]+n[1]+n[2]+n[3]+n[4] == 0) exit
            split(u, p, "/"); printf "%s", p[1]
            for (k = 0; k < 5; k++) printf "\t%.0f\t%d", b[k], n[k]
            printf "\t%.0f\n", newest }'"'"'' > $RAW

# reduce per top-level folder and join with kept.tsv
awk -F'\t' -v now=$NOW -v kept=$A/kept.tsv '
  BEGIN { while ((getline l < kept) > 0) { if (l ~ /^#/ || l ~ /^folder\t/) continue
          split(l, f, "\t"); st[f[1]] = f[2]; nt[f[1]] = f[3] } }
  { t = $1; for (k = 0; k < 5; k++) { B[t,k] += $(2+2*k); N[t,k] += $(3+2*k) }
    if ($12 > NW[t]) NW[t] = $12; seen[t] = 1 }
  END {
    printf "%-30s %-8s %8s %8s %8s %8s %8s %6s  %s\n", "folder", "status", "GB>45d", "<=7d", "8-14d", "15-28d", "safer", "whole", "note"
    printf "%-30s %-8s %8s %8s %8s %8s %8s %6s\n", "", "", "(now)", "", "", "", "", "in(d)"
    for (t in seen) {
      tot = 0; for (k = 0; k < 5; k++) tot += B[t,k]
      w = 45 - (now - NW[t]) / 86400; s = (t in st) ? st[t] : "UNLISTED"
      printf "%-30s %-8s %8.2f %8.2f %8.2f %8.2f %8.2f %6.0f  %s\n", t, s, B[t,0]/1e9, B[t,1]/1e9, B[t,2]/1e9, B[t,3]/1e9, B[t,4]/1e9, w, nt[t]
    } }' $RAW | { read -r h1; read -r h2; echo "$h1"; echo "$h2"; sort -k8,8n; } > $OUT/report_$DAY.body

{
  echo "Scratch purge watch, $(date '+%F %T'). Scratch: $S (45-day purge; age = days since last write)."
  echo "GB>45d: already eligible for deletion. <=7d / 8-14d / 15-28d: becomes eligible within that many days."
  echo "whole in(d): days until even the newest file in the folder is 45 days old (negative = whole folder eligible now)."
  echo "status: KEPT = code in git and/or data archived in $G; REBUILD = deliberately not kept; NOT_KEPT / UNLISTED = decide."
  echo
  echo "ACTION NEEDED (NOT_KEPT or UNLISTED, with data eligible within 14 days):"
  awk 'NR>2 && ($2=="NOT_KEPT" || $2=="UNLISTED") && ($3+$4+$5) > 0 {printf "  %s (%.2f GB within 14 d)\n", $1, $3+$4+$5}' $OUT/report_$DAY.body
  echo
  cat $OUT/report_$DAY.body
} > $OUT/report_$DAY.txt
rm -f $OUT/report_$DAY.body
cp $OUT/report_$DAY.txt $OUT/latest.txt
echo "done $(date '+%F %T'): $OUT/report_$DAY.txt"
