#!/bin/bash
# 47 Tuc isochrones, 8.0-16.0 Gyr in 0.1 Gyr steps, from Tracks_tuc_noco with pyeep (same recipe as audit_ngc6791/build_tuc_iso.sh)
set -e
O=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso_dense
P=/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox/pyeep
T=/dartfs-hpc/scratch/f007058/Tracks_tuc_noco
V=/dartfs-hpc/scratch/f007058/47Tuc_var
RUN=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_fix/dense_run
for R in "$@"; do
  W=$O/w$R; rm -rf $W; mkdir -p $W/raw $W/mesa $W/eeps $O/$R/isos $W/run
  cp $RUN/input.nml $RUN/input.format $RUN/input.eep $RUN/my_history_columns_*.list $W/run/
  tmp=$(mktemp -d); tar xf $T/mcfeh077.$R.tar -C $tmp
  for f in $tmp/m[0-9]*.gz; do [ -s "$f" ] || continue; b=$(basename $f .gz); gunzip -c "$f" > $W/raw/$b; [ -s $W/raw/$b ] || rm -f $W/raw/$b; done
  rm -rf $tmp
  python3 $P/reformat_mesa.py $W/raw -o $W/mesa -j 8 > $W/reformat.log 2>&1
  python3 - "$R" "$W" "$V" "$O" <<'PY'
import os, sys, subprocess
R, W, V, O = sys.argv[1:5]
vals = [l.split()[0] for l in open(f'{V}/varfeh-077.{R}') if l.strip()]
feh, yp, afe = float(vals[0]), float(vals[1]), float(vals[2])
z = float(subprocess.run(['/dartfs-hpc/rc/home/8/f007058/mcdsep/run/gs98getz', str(feh), str(afe)], capture_output=True, text=True).stdout.strip())
tr = sorted(os.listdir(f'{W}/mesa'))
lines = ['#version string, max 8 characters', 'pyeep', '#initial Y, initial Z, [Fe/H], [alpha/Fe], v/vcrit (space separated)',
         '   %.4f %.5e   %.2f        %.2f     0.0' % (yp, z, feh, afe), '#data directories: 1) history files, 2) eeps, 3) isochrones',
         f'{W}/mesa', f'{W}/eeps', f'{O}/{R}/isos', '#read history_columns', 'my_history_columns_basic.list',
         '#specify tracks', str(len(tr))] + tr + ['#specify isochrones', 'isochrones_example.txt', 'min_max', 'linear', '81', '8.0e9', '16.0e9']
open(f'{W}/run/in_tuc{R}', 'w').write('\n'.join(lines) + '\n')
print(f'  {R}: Y={yp:.4f} Z={z:.5e} [Fe/H]={feh:.3f} [a/Fe]={afe:.2f}, {len(tr)} tracks')
PY
  cd $W/run
  python3 $P/make_eep.py in_tuc$R -j 4 > $W/eep.log 2>&1
  python3 $P/make_iso.py in_tuc$R > $W/iso.log 2>&1
  echo "  $R -> $(stat -c %s $O/$R/isos/isochrones_example.txt) bytes"
done
# lower-MS patch, as build_tuc_one.sh
for R in "$@"; do
  mv /dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso_dense/$R/isos/isochrones_example.txt /dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso_dense/$R/isos/isochrones_pyeep.txt
  /dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python /dartfs-hpc/scratch/f007058/scmd_tb/patch_lower_ms.py /dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso_dense/$R/isos/isochrones_pyeep.txt /dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso_dense/w$R/mesa /dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso_dense/$R/isos/isochrones_example.txt
done
