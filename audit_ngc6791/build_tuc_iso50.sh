#!/bin/bash
A=/dartfs-hpc/scratch/f007058/audit_ngc6791
W=$A/pyeep_work
P=/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox-main/pyeep
T=/dartfs-hpc/scratch/f007058/Tracks_tuc_noco
V=/dartfs-hpc/scratch/f007058/47Tuc_var
for R in $(cat $A/tuc_good50.txt); do
  [ -f $W/isos/tuc$R/isochrones_tuc$R.txt ] && { echo "$R done"; continue; }
  rm -rf $W/raw/tuc$R $W/mesa/tuc$R; mkdir -p $W/raw/tuc$R $W/mesa/tuc$R $W/eeps/tuc$R $W/isos/tuc$R
  rm -rf $A/extracted/tuc$R; mkdir -p $A/extracted/tuc$R
  tar xf $T/mcfeh077.$R.tar -C $A/extracted/tuc$R
  for f in $A/extracted/tuc$R/m[0-9]*.gz; do
    [ -s "$f" ] || continue
    b=$(basename $f .gz); gunzip -c "$f" > $W/raw/tuc$R/$b
    [ -s $W/raw/tuc$R/$b ] || rm -f $W/raw/tuc$R/$b
  done
  python3 $P/reformat_mesa.py $W/raw/tuc$R -o $W/mesa/tuc$R -j 8 >/dev/null 2>&1
  python3 - <<PY
import os,subprocess
W="$W"; R="$R"
vals=[l.split()[0] for l in open("$V/varfeh-077.$R") if l.strip()]
feh,yp,afe=float(vals[0]),float(vals[1]),float(vals[2])
z=float(subprocess.run(['/dartfs-hpc/rc/home/8/f007058/mcdsep/run/gs98getz',str(feh),str(afe)],
                       capture_output=True,text=True).stdout.strip())
tr=sorted(os.listdir(f'{W}/mesa/tuc{R}'))
lines=['#version string, max 8 characters','pyeep',
 '#initial Y, initial Z, [Fe/H], [alpha/Fe], v/vcrit (space separated)',
 '   %.4f %.5e   %.2f        %.2f     0.0'%(yp,z,feh,afe),
 '#data directories: 1) history files, 2) eeps, 3) isochrones',
 f'{W}/mesa/tuc{R}', f'{W}/eeps/tuc{R}', f'{W}/isos/tuc{R}',
 '#read history_columns','my_history_columns_basic.list',
 '#specify tracks',str(len(tr))]+tr+[
 '#specify isochrones',f'isochrones_tuc{R}.txt','min_max','linear','14','6.0e9','8.6e9']
open(f'{W}/run/in_tuc{R}','w').write('\n'.join(lines)+'\n')
print('  %s: Y=%.4f Z=%.5e feh=%.3f afe=%.2f tracks=%d'%(R,yp,z,feh,afe,len(tr)))
PY
  cd $W/run
  python3 $P/make_eep.py in_tuc$R -j 8 > eep_tuc$R.log 2>&1
  python3 $P/make_iso.py in_tuc$R > iso_tuc$R.log 2>&1
  echo "    -> $(ls -la $W/isos/tuc$R/isochrones_tuc$R.txt 2>/dev/null | awk '{print $5}') bytes"
  rm -rf $W/raw/tuc$R
done
echo "DONE: $(ls $W/isos/tuc*/isochrones_tuc*.txt 2>/dev/null | wc -l) isochrone files"
