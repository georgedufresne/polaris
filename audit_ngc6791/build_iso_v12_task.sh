#!/bin/bash
# Build isochrones under the pyeep v1.2 DSEE config (input_dsee.nml + input_dsee.eep).
# Writes to *_v12 directories so the existing v1.0 isochrones are left intact.
# $1 = chunk index (0-based), $2 = chunk size
A=/dartfs-hpc/scratch/f007058/audit_ngc6791
W=$A/pyeep_work
P=/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox/pyeep
T=/dartfs-hpc/scratch/f007058/Tracks_noco
V=/dartfs-hpc/rc/home/8/f007058/mcdsep/MC_var_noco
CHUNK=$1; SIZE=$2
RD=$W/run_v12_t$CHUNK; export RD; rm -rf $RD; mkdir -p $RD
# v1.2: stage input_dsee.nml AS input.nml (it references input_dsee.eep by name)
cp $P/run/input_dsee.nml $RD/input.nml
cp $P/run/input_dsee.eep $RD/input_dsee.eep
cp $P/run/input.format $P/run/my_history_columns_*.list $RD/
start=$(( CHUNK*SIZE + 1 )); end=$(( start + SIZE - 1 ))
for R in $(sed -n "${start},${end}p" $A/v12_runs.txt); do
  [ -s $W/isos_v12/v12$R/isochrones_v12$R.txt ] && { echo "$R exists"; continue; }
  mkdir -p $W/raw_v12_t$CHUNK $W/mesa_v12/v12$R $W/eeps_v12/v12$R $W/isos_v12/v12$R
  rm -rf $W/raw_v12_t$CHUNK/*; mkdir -p $W/raw_v12_t$CHUNK/$R
  tmp=$(mktemp -d); tar xf $T/mcfeh031.$R.tar -C $tmp 2>/dev/null
  for f in $tmp/m[0-9]*.gz; do
    [ -s "$f" ] || continue
    b=$(basename $f .gz); gunzip -c "$f" > $W/raw_v12_t$CHUNK/$R/$b
    [ -s $W/raw_v12_t$CHUNK/$R/$b ] || rm -f $W/raw_v12_t$CHUNK/$R/$b
  done
  rm -rf $tmp
  python3 $P/reformat_mesa.py $W/raw_v12_t$CHUNK/$R -o $W/mesa_v12/v12$R -j 1 >/dev/null 2>&1
  python3 - "$R" "$W" "$V" <<'PY'
import os,sys,subprocess
R,W,V=sys.argv[1],sys.argv[2],sys.argv[3]
vals=[l.split()[0] for l in open(f'{V}/varfeh+031.{R}') if l.strip()]
feh,yp,afe=float(vals[0]),float(vals[1]),float(vals[2])
z=float(subprocess.run(['/dartfs-hpc/rc/home/8/f007058/mcdsep/run/gs98getz',str(feh),str(afe)],
                       capture_output=True,text=True).stdout.strip())
tr=sorted(os.listdir(f'{W}/mesa_v12/v12{R}'))
lines=['#version string, max 8 characters','pyeep',
 '#initial Y, initial Z, [Fe/H], [alpha/Fe], v/vcrit (space separated)',
 '   %.4f %.5e   %.2f        %.2f     0.0'%(yp,z,feh,afe),
 '#data directories: 1) history files, 2) eeps, 3) isochrones',
 f'{W}/mesa_v12/v12{R}', f'{W}/eeps_v12/v12{R}', f'{W}/isos_v12/v12{R}',
 '#read history_columns','my_history_columns_basic.list',
 '#specify tracks',str(len(tr))]+tr+[
 '#specify isochrones',f'isochrones_v12{R}.txt','min_max','linear','81','4.0e9','12.0e9']
open(os.environ['RD']+f'/in_v12{R}','w').write('\n'.join(lines)+'\n')
PY
  cd $RD
  python3 $P/make_eep.py in_v12$R -j 1 > $W/logs/eep_v12$R.log 2>&1
  python3 $P/make_iso.py in_v12$R > $W/logs/iso_v12$R.log 2>&1
  sz=$(stat -c %s $W/isos_v12/v12$R/isochrones_v12$R.txt 2>/dev/null || echo 0)
  echo "$R -> $sz bytes"
  rm -rf $W/raw_v12_t$CHUNK/$R $W/mesa_v12/v12$R
done
rm -rf $W/raw_v12_t$CHUNK $RD
