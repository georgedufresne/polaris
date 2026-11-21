#!/bin/bash
# Tar the per-mass dsep4 tracks into per-run archives, then build v1.2 isochrones.
# $1 = run number
set -u
R=$1
VER=${2:-v12}
PARTS=/dartfs-hpc/scratch/f007058/Tracks_dsep4_parts/$R
STORE=/dartfs-hpc/scratch/f007058/Tracks_dsep4
P=/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox/pyeep
W=/dartfs-hpc/scratch/f007058/dsep4_batch/pyeep_work
V=/dartfs-hpc/rc/home/8/f007058/mcdsep/MC_var_noco_d4
mkdir -p $STORE $W/logs $W/mesa/${VER}_d4$R $W/eeps/${VER}_d4$R $W/isos/${VER}_d4$R

n=$(ls $PARTS/m[0-9]*.gz 2>/dev/null | wc -l)
echo "run $R: $n mass tracks present"
[ "$n" -eq 0 ] && { echo "  nothing to do"; exit 1; }

( cd $PARTS && tar -cf $STORE/mcfeh031.$R.tar m[0-9]*.gz \
                && tar -cf $STORE/iso_mcfeh031.$R.tar iso_m[0-9]*.gz )

# unpack tracks for pyeep
RAW=$W/raw/${VER}_$R; rm -rf $RAW; mkdir -p $RAW
for f in $PARTS/m[0-9]*.gz; do
  b=$(basename $f .gz); gunzip -c "$f" > $RAW/$b
  [ -s $RAW/$b ] || rm -f $RAW/$b
done
# dsep4 writes native MESA (row-per-model) history directly, unlike dsepX's
# one-value-per-line output, so reformat_mesa.py correctly skips every file and
# writes nothing. Convert what needs it, then copy across anything it skipped.
python3 $P/reformat_mesa.py $RAW -o $W/mesa/${VER}_d4$R -j 4 > $W/logs/reformat_$R.log 2>&1
for f in $RAW/*; do
  b=$(basename "$f")
  [ -s "$W/mesa/${VER}_d4$R/$b" ] || cp -f "$f" "$W/mesa/${VER}_d4$R/$b"
done
echo "  mesa dir: $(ls $W/mesa/${VER}_d4$R | wc -l) tracks"

# v1.2 DSEE config: input_dsee.nml staged AS input.nml (it names input_dsee.eep)
RD=$W/run_${VER}_$R; rm -rf $RD; mkdir -p $RD; export RD
if [ "$VER" = "v10" ]; then
  cp $P/run/input.nml $RD/input.nml
  cp $P/run/input.eep $RD/input.eep
else
  cp $P/run/input_dsee.nml $RD/input.nml
  cp $P/run/input_dsee.eep $RD/input_dsee.eep
fi
cp $P/run/input.format $P/run/my_history_columns_*.list $RD/

python3 - "$R" "$W" "$V" "$VER" <<'PY'
import os,sys,subprocess
R,W,V2,V=sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4]
vals=[l.split()[0] for l in open(f'{V2}/varfeh+031.{R}') if l.strip()]
feh,yp,afe=float(vals[0]),float(vals[1]),float(vals[2])
z=float(subprocess.run(['/dartfs-hpc/rc/home/8/f007058/mcdsep/run/gs98getz',str(feh),str(afe)],
                       capture_output=True,text=True).stdout.strip())
tr=sorted(os.listdir(f'{W}/mesa/{V}_d4{R}'))
lines=['#version string, max 8 characters','pyeep',
 '#initial Y, initial Z, [Fe/H], [alpha/Fe], v/vcrit (space separated)',
 '   %.4f %.5e   %.2f        %.2f     0.0'%(yp,z,feh,afe),
 '#data directories: 1) history files, 2) eeps, 3) isochrones',
 f'{W}/mesa/{V}_d4{R}', f'{W}/eeps/{V}_d4{R}', f'{W}/isos/{V}_d4{R}',
 '#read history_columns','my_history_columns_basic.list',
 '#specify tracks',str(len(tr))]+tr+[
 '#specify isochrones',f'isochrones_{V}d4{R}.txt','min_max','linear','81','4.0e9','12.0e9']
open(os.environ['RD']+f'/in_d4{R}','w').write('\n'.join(lines)+'\n')
PY

cd $RD
python3 $P/make_eep.py in_${VER}d4$R -j 4 > $W/logs/eep_${VER}d4$R.log 2>&1
python3 $P/make_iso.py in_${VER}d4$R      > $W/logs/iso_${VER}d4$R.log 2>&1
sz=$(stat -c %s $W/isos/${VER}_d4$R/isochrones_${VER}d4$R.txt 2>/dev/null || echo 0)
echo "run $R -> isochrone $sz bytes"
rm -rf $RAW $W/mesa/${VER}_d4$R $RD
