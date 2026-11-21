#!/bin/bash
A=/dartfs-hpc/scratch/f007058/audit_ngc6791
W=$A/pyeep_work
P=/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox-main/pyeep
E=$A/extracted
for V in $(cat /scratch/f007058/claude-550484/-dartfs-hpc/84db2d29-6b4d-4b69-a72e-63390c149db0/scratchpad/good_vars.txt); do
  [ -f $W/isos/$V/isochrones_$V.txt ] && { echo "$V already done"; continue; }
  rm -rf $W/raw/$V $W/mesa/$V; mkdir -p $W/raw/$V $W/mesa/$V
  for f in $E/$V/m[0-9]*.gz; do b=$(basename $f .gz); gunzip -c "$f" > $W/raw/$V/$b; done
  python3 $P/reformat_mesa.py $W/raw/$V -o $W/mesa/$V -j 8 >/dev/null 2>&1
  python3 $A/mk_input.py $V
  cd $W/run
  python3 $P/make_eep.py in_$V -j 8 >$W/run/eep_$V.log 2>&1
  python3 $P/make_iso.py in_$V >$W/run/iso_$V.log 2>&1
  echo "  $V -> $(ls -la $W/isos/$V/isochrones_$V.txt 2>/dev/null | awk '{print $5}') bytes"
  rm -rf $W/raw/$V
done
echo "ALL DONE: $(ls $W/isos/*/isochrones_*.txt 2>/dev/null | wc -l) isochrone files"
