#!/bin/bash
# Build isochrones for a chunk of runs. $1 = chunk index (0-based), $2 = chunk size
A=/dartfs-hpc/scratch/f007058/audit_ngc6791
W=$A/pyeep_work
P=/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox-main/pyeep
T=/dartfs-hpc/scratch/f007058/Tracks_tuc_noco
V=/dartfs-hpc/rc/home/8/f007058/mcdsep/MC_var_tuc
CHUNK=$1; SIZE=$2
# private run dir so concurrent make_eep calls don't collide on control files
RD=$W/run_t$CHUNK; export RD; rm -rf $RD; mkdir -p $RD
cp $W/run/input.nml $W/run/input.format $W/run/input.eep $W/run/my_history_columns_*.list $RD/
start=$(( CHUNK*SIZE + 1 )); end=$(( start + SIZE - 1 ))
for R in $(sed -n "${start},${end}p" $A/tuc_all.txt); do
  [ -s $W/isos/tucOLD$R/isochrones_tucOLD$R.txt ] && { echo "$R exists"; continue; }
  mkdir -p $W/raw_t$CHUNK $W/mesa/tuc$R $W/eeps/tucOLD$R $W/isos/tucOLD$R
  rm -rf $W/raw_t$CHUNK/*; mkdir -p $W/raw_t$CHUNK/$R
  tmp=$(mktemp -d); tar xf $T/mcfeh077.$R.tar -C $tmp 2>/dev/null
  for f in $tmp/m[0-9]*.gz; do
    [ -s "$f" ] || continue
    b=$(basename $f .gz); gunzip -c "$f" > $W/raw_t$CHUNK/$R/$b
    [ -s $W/raw_t$CHUNK/$R/$b ] || rm -f $W/raw_t$CHUNK/$R/$b
  done
  rm -rf $tmp
  python3 $P/reformat_mesa.py $W/raw_t$CHUNK/$R -o $W/mesa/tuc$R -j 1 >/dev/null 2>&1
  python3 - "$R" "$W" "$V" <<'PY'
import os,sys,subprocess
R,W,V=sys.argv[1],sys.argv[2],sys.argv[3]
vals=[l.split()[0] for l in open(f'{V}/varfeh-077.{R}') if l.strip()]
feh,yp,afe=float(vals[0]),float(vals[1]),float(vals[2])
z=float(subprocess.run(['/dartfs-hpc/rc/home/8/f007058/mcdsep/run/gs98getz',str(feh),str(afe)],
                       capture_output=True,text=True).stdout.strip())
tr=sorted(os.listdir(f'{W}/mesa/tuc{R}'))
lines=['#version string, max 8 characters','pyeep',
 '#initial Y, initial Z, [Fe/H], [alpha/Fe], v/vcrit (space separated)',
 '   %.4f %.5e   %.2f        %.2f     0.0'%(yp,z,feh,afe),
 '#data directories: 1) history files, 2) eeps, 3) isochrones',
 f'{W}/mesa/tuc{R}', f'{W}/eeps/tucOLD{R}', f'{W}/isos/tucOLD{R}',
 '#read history_columns','my_history_columns_basic.list',
 '#specify tracks',str(len(tr))]+tr+[
 '#specify isochrones',f'isochrones_tucOLD{R}.txt','min_max','linear','26','10.0e9','15.0e9']
open(os.environ['RD']+f'/in_tuc{R}','w').write('\n'.join(lines)+'\n')
PY
  cd $RD
  python3 $P/make_eep.py in_tuc$R -j 1 > /dev/null 2>&1
  python3 $P/make_iso.py in_tuc$R > /dev/null 2>&1
  sz=$(stat -c %s $W/isos/tucOLD$R/isochrones_tucOLD$R.txt 2>/dev/null || echo 0)
  echo "$R -> $sz bytes"
  rm -rf $W/raw_t$CHUNK/$R $W/mesa/tuc$R
done
rm -rf $W/raw_t$CHUNK $RD
