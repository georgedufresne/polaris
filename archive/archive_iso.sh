#!/bin/bash
#SBATCH --job-name=archive_iso
#SBATCH --account=free
#SBATCH --partition=standard
#SBATCH --cpus-per-task=8
#SBATCH --ntasks=1
#SBATCH --time=12:00:00
#SBATCH --output=/dartfs-hpc/rc/lab/C/ChaboyerB/george/polaris/archive/logs/archive_iso_%j.log
# Copy the physics-on (Monte Carlo) model outputs from Polaris scratch to the lab share before the 45-day purge.
#   iso/        DSEP per-mass iso files (iso_m*.feh*.gz), one tar per run, as DSEP3 writes them (rebuild isochrones with pyeep)
#   isochrones/ the built isochrones the sCMD pipeline reads (<run>/isos/isochrones_example.txt), one .tar.xz per set
#   var/        the Monte Carlo var files of each grid
#   obs/        observational inputs the fits read (M15 and 47 Tuc star catalogs and artificial-star tests)
# Usage: sbatch archive_iso.sh [item ...]   (default: every item except tuc4, which waits for the 47 Tuc DSEP4 grid)
# Copies only; never deletes from scratch. Each item is checked (file counts, archive integrity) and logged.
set -uo pipefail
S=/dartfs-hpc/scratch/f007058
G=/dartfs-hpc/rc/lab/C/ChaboyerB/george
LOG=$G/polaris/archive/archive_log.tsv
MINFREE_GB=20
T=${SLURM_CPUS_PER_TASK:-4}

say()  { echo "$(date '+%F %T') $*"; }
free_ok() {
  local f; f=$(df -BG --output=avail $G | tail -1 | tr -dc 0-9)
  [ "$f" -ge $MINFREE_GB ] || { say "STOP: lab share has ${f} GB free (< $MINFREE_GB)"; exit 1; }
}
record() { echo -e "$(date '+%F %T')\t$1\t$2\t$3\t$4" >> $LOG; say "$1: $3 -> $2 ($4)"; }

# copy DSEP3 per-run iso tars as they are (contents already gzipped)
copy_tars() {   # name src_dir glob
  local name=$1 src=$2 glob=$3 dst=$G/iso/$1; free_ok; mkdir -p $dst
  rsync -rlt --include="$glob" --exclude='*' $src/ $dst/
  local ns nd bad=0
  ns=$(find $src -maxdepth 1 -name "$glob" | wc -l); nd=$(find $dst -maxdepth 1 -name "$glob" | wc -l)
  for f in $dst/$glob; do tar -tf $f > /dev/null 2>&1 || { bad=$((bad+1)); say "bad tar $f"; }; done
  record $name $dst "$src/$glob" "$nd of $ns tars, $bad unreadable, $(du -sh $dst | cut -f1)"
}

# pack DSEP4 per-run iso files (run dirs of iso_m*.gz) into one tar per run, same naming as DSEP3
pack_runs() {   # name src_dir tag
  local name=$1 src=$2 tag=$3 dst=$G/iso/$1 n=0 bad=0; free_ok; mkdir -p $dst
  for d in $src/[0-9]*/; do
    r=$(basename $d); ls $d/iso_m*.gz > /dev/null 2>&1 || continue
    out=$dst/iso_mc$tag.$r.tar
    if [ ! -s $out ]; then (cd $d && tar -cf $out.tmp iso_m*.gz) && mv $out.tmp $out; fi
    a=$(tar -tf $out | wc -l); b=$(ls $d/iso_m*.gz | wc -l)
    [ "$a" -eq "$b" ] || { bad=$((bad+1)); say "count mismatch $r: $a vs $b"; }
    n=$((n+1))
  done
  record $name $dst "$src/<run>/iso_m*.gz" "$n runs packed, $bad mismatched, $(du -sh $dst | cut -f1)"
}

# built isochrones: one .tar.xz of <set>/<run>/isos/isochrones_example.txt
pack_isochrones() {   # name parent set_dir [pattern]
  local name=$1 parent=$2 set=$3 pat=${4:-isochrones_example.txt} out=$G/isochrones/$1.tar.xz; free_ok
  local ns na
  # -L / -h: follow symlinks and store the file contents (Isochrones_noco holds links into audit_ngc6791/pyeep_work)
  ns=$(cd $parent && find -L $set -name "$pat" -type f | wc -l)
  (cd $parent && find -L $set -name "$pat" -type f -print0 | sort -z | tar -h --null -T - -cf -) | xz -T$T -6 > $out.tmp \
    && mv $out.tmp $out
  na=$(xz -T$T -dc $out | tar -tf - | wc -l)
  # spot check: first file in the archive is byte-identical to the source
  f1=$(xz -dc $out | tar -tf - | head -1)
  same=$(cmp -s <(xz -dc $out | tar -xOf - "$f1") $parent/$f1 && echo identical || echo DIFFERENT)
  record $name $out "$parent/$set/.../$pat" "$na of $ns files, spot check $f1 $same, $(du -sh $out | cut -f1)"
}

# observational inputs (catalogs, artificial-star tests): one .tar.xz of the listed paths (relative to parent), as they are
pack_paths() {   # name parent path ...
  local name=$1 parent=$2; shift 2; local out=$G/obs/$name.tar.xz ns na f1 same; free_ok; mkdir -p $G/obs
  ns=$(cd $parent && find "$@" -type f | wc -l)
  (cd $parent && find "$@" -type f -print0 | sort -z | tar --null -T - -cf -) | xz -T$T -6 > $out.tmp && mv $out.tmp $out
  na=$(xz -T$T -dc $out | tar -tf - | wc -l)
  f1=$(xz -dc $out | tar -tf - | head -1)
  same=$(cmp -s <(xz -dc $out | tar -xOf - "$f1") $parent/$f1 && echo identical || echo DIFFERENT)
  record $name $out "$parent: $*" "$na of $ns files, spot check $f1 $same, $(du -sh $out | cut -f1)"
}

pack_var() {   # name src_dir
  local out=$G/var/$1.tar.xz; free_ok
  (cd $(dirname $2) && tar -cf - $(basename $2)) | xz -T$T -6 > $out.tmp && mv $out.tmp $out
  record $1 $out "$2" "$(xz -dc $out | tar -tf - | grep -vc '/$') files, $(du -sh $out | cut -f1)"
}

ITEMS=${*:-var ngc6791 tuc3 m15 tuc3_10k}
mkdir -p $G/polaris/archive/logs
[ -f $LOG ] || echo -e "time\titem\tdestination\tsource\tcheck" > $LOG
say "items: $ITEMS"
for it in $ITEMS; do
  case $it in
    var)      pack_var var_47tuc $S/47Tuc_var; pack_var var_m15 $S/M15/dsep4/var; pack_var var_ngc6791 $S/emul_W14/var
              cp -p $S/M15/M15_var.tar $G/var/ ;;
    ngc6791)  copy_tars ngc6791_dsep3_mc $S/Tracks_noco 'iso_mcfeh031.*.tar' ;;
    ngc6791iso) pack_isochrones ngc6791_dsep3_mc $S Isochrones_noco ;;
    tuc3)     copy_tars 47tuc_dsep3_mc $S/Tracks_tuc_noco 'iso_mcfeh077.*.tar'
              pack_isochrones 47tuc_dsep3_mc $S/scmd_tb tuc_iso ;;
    m15)      pack_runs m15_dsep4_mc $S/M15/tracks_dsep4mc feh227
              pack_isochrones m15_dsep4_mc $S/M15 iso_dsep4mc ;;
    tuc3_10k) pack_isochrones 47tuc_dsep3_10k_dsepformat $S 47Tuc_iso 'feh077*' ;;
    obs)      pack_paths obs_m15_47tuc $S Files_from_Martin/47Tuc_data Files_from_Martin/M15_data Files_from_Martin/NGC0104.47T \
                Files_from_Martin/Photometry/47Tuc_inputfiles Files_from_Martin/Photometry/M15_inputfiles \
                M15/NGC7078.M15 scmd_tb/m15_grid_d4/inputs scmd_tb/hugs/47Tuc_hugs_meth2.dat ;;
    tuc4)     pack_runs 47tuc_dsep4_mc $S/Tracks_tuc_dsep4_mcfix feh077
              pack_isochrones 47tuc_dsep4_mc $S/scmd_tb tuc_iso_d4v2_mcfix ;;
    *)        say "unknown item $it" ;;
  esac
done
say "done; lab share free: $(df -h --output=avail $G | tail -1)"
