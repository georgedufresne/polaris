#!/bin/bash
# One M15 DSEP4 (Monte Carlo ON) run -> pyeep isochrones 8.0-16.0 Gyr in 0.1 Gyr steps,
# main-sequence holes filled from the tracks. As tuc_dsep4_batch/build_tuc_d4_iso.sh, except:
#  - pyeep's iso_controls stop the isochrones at EEP 605 (RGB tip): set_max_eep_number
#    = .true., new_max_eep_number = 605 (nml/input.nml). pyeep drops every low-mass track
#    above 0.6 Msun with fewer EEPs than the longest one; DSEP4 takes some masses past the
#    He flash (EEP 631) and stops others at the tip (605), so the 605 tracks were dropped
#    and isochrones older than ~10.5 Gyr lost their turnoff. No DSEP3 track passes 605,
#    so for DSEP3 this setting changes nothing and the two sets stay comparable;
#  - everything it reads is on scratch (no Kerberos ticket needed): pyeep, the var
#    files, gs98getz and the namelists are copies in this directory;
#  - M15 copy (2026-11-02): tracks M15/tracks_dsep4, var files M15/dsep4/var, output
#    M15/iso_dsep4; pyeep, its namelists and gs98getz are the 47 Tuc build copies.
#   $1 = run number. T, O may be overridden from the environment.
set -u
R=$1
D=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_d4v2_build
T=${T:-/dartfs-hpc/scratch/f007058/M15/tracks_dsep4}
O=${O:-/dartfs-hpc/scratch/f007058/M15/iso_dsep4}
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
P=$D/pyeep
PATCH=/dartfs-hpc/scratch/f007058/scmd_tb/patch_lower_ms.py
F=$O/$R/isos/isochrones_example.txt
[ -s $O/$R/isos/isochrones_pyeep.txt ] && [ -s $F ] && { echo "$R iso done"; exit 0; }
n=$(ls $T/$R/m[0-9]*.gz 2>/dev/null | wc -l)
echo "run $R: $n mass tracks"
[ "$n" -eq 0 ] && { echo "  no tracks"; exit 1; }
( cd $T/$R && tar -cf $T/mcfeh227.$R.tar m[0-9]*.gz && tar -cf $T/iso_mcfeh227.$R.tar iso_m[0-9]*.gz )
W=$O/w$R; rm -rf "${W:?}"; mkdir -p $W/raw $W/mesa $W/eeps $O/$R/isos $W/run
cp $D/nml/input.nml $D/nml/input.format $D/nml/input.eep $D/nml/my_history_columns_*.list $W/run/
for f in $T/$R/m[0-9]*.gz; do
  b=$(basename $f .gz); gunzip -c "$f" > $W/raw/$b
  [ -s $W/raw/$b ] || rm -f $W/raw/$b
done
$PY $P/reformat_mesa.py $W/raw -o $W/mesa -j 1 > $W/reformat.log 2>&1
for f in $W/raw/*; do b=$(basename "$f"); [ -s "$W/mesa/$b" ] || cp -f "$f" "$W/mesa/$b"; done
$PY - "$R" "$W" "/dartfs-hpc/scratch/f007058/M15/dsep4/var/varfeh-227.$R" "$O" "$D/gs98getz" <<'PY'
import os, sys, subprocess
R, W, VF, O, GETZ = sys.argv[1:6]
vals = [l.split()[0] for l in open(VF) if l.strip()]
feh, yp, afe = float(vals[0]), float(vals[1]), float(vals[2])
z = float(subprocess.run([GETZ, str(feh), str(afe)], capture_output=True, text=True).stdout.strip())
tr = sorted(os.listdir(f'{W}/mesa'))
lines = ['#version string, max 8 characters', 'pyeep', '#initial Y, initial Z, [Fe/H], [alpha/Fe], v/vcrit (space separated)',
         '   %.4f %.5e   %.2f        %.2f     0.0' % (yp, z, feh, afe), '#data directories: 1) history files, 2) eeps, 3) isochrones',
         f'{W}/mesa', f'{W}/eeps', f'{O}/{R}/isos', '#read history_columns', 'my_history_columns_basic.list',
         '#specify tracks', str(len(tr))] + tr + ['#specify isochrones', 'isochrones_example.txt', 'min_max', 'linear', '81', '8.0e9', '16.0e9']
open(f'{W}/run/in_tuc{R}', 'w').write('\n'.join(lines) + '\n')
print(f'  {R}: Y={yp:.4f} Z={z:.5e} [Fe/H]={feh:.3f} [a/Fe]={afe:.2f}, {len(tr)} tracks')
PY
cd $W/run
$PY $P/make_eep.py in_tuc$R -j 1 > $W/eep.log 2>&1
$PY $P/make_iso.py in_tuc$R > $W/iso.log 2>&1
[ -s $F ] || { echo "  $R ISOCHRONE BUILD FAILED (see $W/eep.log, $W/iso.log)"; exit 1; }
mv $F $O/$R/isos/isochrones_pyeep.txt
$PY $PATCH $O/$R/isos/isochrones_pyeep.txt $W/mesa $F
rm -rf "${W:?}/raw"
echo "  $R -> $(stat -c %s $F) bytes"
