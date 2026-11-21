#!/bin/bash
# One 47 Tuc DSEP4 (Monte Carlo ON) run: tar the per-mass tracks, then build pyeep
# isochrones 8.0-16.0 Gyr in 0.1 Gyr steps and fill main-sequence holes from the
# tracks. Same recipe as scmd_tb/build_tuc_12gyr.sh + build_tuc_one.sh (the DSEP3
# 47 Tuc set in scmd_tb/tuc_iso), so the two sets compare one-to-one.
#   $1 = run number. T, V, O, VARF may be overridden from the environment.
set -u
R=$1
T=${T:-/dartfs-hpc/scratch/f007058/Tracks_tuc_dsep4_mc}
V=${V:-/dartfs-hpc/rc/home/8/f007058/mcdsep/MC_var_tuc_d4}
VARF=${VARF:-varfeh-077.$R}
O=${O:-/dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso_d4}
P=/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox/pyeep
RUN=/dartfs-hpc/scratch/f007058/audit_ngc6791/pyeep_work/run
PATCH=/dartfs-hpc/scratch/f007058/scmd_tb/patch_lower_ms.py
F=$O/$R/isos/isochrones_example.txt

[ -s $O/$R/isos/isochrones_pyeep.txt ] && [ -s $F ] && { echo "$R iso done"; exit 0; }
n=$(ls $T/$R/m[0-9]*.gz 2>/dev/null | wc -l)
echo "run $R: $n mass tracks"
[ "$n" -eq 0 ] && { echo "  no tracks"; exit 1; }

# per-run archives, same layout as Tracks_tuc_noco
( cd $T/$R && tar -cf $T/mcfeh077.$R.tar m[0-9]*.gz && tar -cf $T/iso_mcfeh077.$R.tar iso_m[0-9]*.gz )

W=$O/w$R; rm -rf $W; mkdir -p $W/raw $W/mesa $W/eeps $O/$R/isos $W/run
cp $RUN/input.nml $RUN/input.format $RUN/input.eep $RUN/my_history_columns_*.list $W/run/
for f in $T/$R/m[0-9]*.gz; do
  b=$(basename $f .gz); gunzip -c "$f" > $W/raw/$b
  [ -s $W/raw/$b ] || rm -f $W/raw/$b
done
# dsep4 already writes MESA row-per-model histories, so reformat_mesa.py skips
# them; copy across whatever it did not write.
python3 $P/reformat_mesa.py $W/raw -o $W/mesa -j 1 > $W/reformat.log 2>&1
for f in $W/raw/*; do b=$(basename "$f"); [ -s "$W/mesa/$b" ] || cp -f "$f" "$W/mesa/$b"; done

python3 - "$R" "$W" "$V/$VARF" "$O" <<'PY'
import os, sys, subprocess
R, W, VF, O = sys.argv[1:5]
vals = [l.split()[0] for l in open(VF) if l.strip()]
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
python3 $P/make_eep.py in_tuc$R -j 1 > $W/eep.log 2>&1
python3 $P/make_iso.py in_tuc$R > $W/iso.log 2>&1
[ -s $F ] || { echo "  $R ISOCHRONE BUILD FAILED (see $W/eep.log, $W/iso.log)"; exit 1; }
mv $F $O/$R/isos/isochrones_pyeep.txt
python3 $PATCH $O/$R/isos/isochrones_pyeep.txt $W/mesa $F
rm -rf $W/raw
echo "  $R -> $(stat -c %s $F) bytes"
