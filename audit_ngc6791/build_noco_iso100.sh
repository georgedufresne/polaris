#!/bin/bash
A=/dartfs-hpc/scratch/f007058/audit_ngc6791
W=$A/pyeep_work
P=/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox-main/pyeep
T=/dartfs-hpc/scratch/f007058/Tracks_noco
V=/dartfs-hpc/rc/home/8/f007058/mcdsep/MC_var_noco
for R in $(cat /dartfs-hpc/scratch/f007058/audit_ngc6791/noco100.txt); do
  [ -f $W/isos/noco$R/isochrones_noco$R.txt ] && { echo "$R done"; continue; }
  rm -rf $W/raw/noco$R $W/mesa/noco$R; mkdir -p $W/raw/noco$R $W/mesa/noco$R $W/eeps/noco$R $W/isos/noco$R
  rm -rf $A/extracted/noco$R; mkdir -p $A/extracted/noco$R
  tar xf $T/mcfeh031.$R.tar -C $A/extracted/noco$R
  for f in $A/extracted/noco$R/m[0-9]*.gz; do
    [ -s "$f" ] || continue
    b=$(basename $f .gz); gunzip -c "$f" > $W/raw/noco$R/$b
    [ -s $W/raw/noco$R/$b ] || rm -f $W/raw/noco$R/$b
  done
  python3 $P/reformat_mesa.py $W/raw/noco$R -o $W/mesa/noco$R -j 8 >/dev/null 2>&1
  python3 - <<PY
import os,subprocess
W="$W"; R="$R"
vals=[l.split()[0] for l in open("$V/varfeh+031.$R") if l.strip()]
feh,yp,afe=float(vals[0]),float(vals[1]),float(vals[2])
z=float(subprocess.run(['/dartfs-hpc/rc/home/8/f007058/mcdsep/run/gs98getz',str(feh),str(afe)],
                       capture_output=True,text=True).stdout.strip())
tr=sorted(os.listdir(f'{W}/mesa/noco{R}'))
lines=['#version string, max 8 characters','pyeep',
 '#initial Y, initial Z, [Fe/H], [alpha/Fe], v/vcrit (space separated)',
 '   %.4f %.5e   %.2f        %.2f     0.0'%(yp,z,feh,afe),
 '#data directories: 1) history files, 2) eeps, 3) isochrones',
 f'{W}/mesa/noco{R}', f'{W}/eeps/noco{R}', f'{W}/isos/noco{R}',
 '#read history_columns','my_history_columns_basic.list',
 '#specify tracks',str(len(tr))]+tr+[
 '#specify isochrones',f'isochrones_noco{R}.txt','min_max','linear','81','4.0e9','12.0e9']
open(f'{W}/run/in_noco{R}','w').write('\n'.join(lines)+'\n')
print('  %s: tracks=%d'%(R,len(tr)))
PY
  cd $W/run
  python3 $P/make_eep.py in_noco$R -j 8 > eep_noco$R.log 2>&1
  python3 $P/make_iso.py in_noco$R > iso_noco$R.log 2>&1
  echo "    -> $(ls -la $W/isos/noco$R/isochrones_noco$R.txt 2>/dev/null | awk '{print $5}') bytes"
  rm -rf $W/raw/noco$R
done
echo "DONE: $(ls $W/isos/noco*/isochrones_noco*.txt 2>/dev/null | wc -l)"
