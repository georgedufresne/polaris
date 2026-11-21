#!/bin/bash
# one 47 Tuc run: pyeep isochrones 8-16 Gyr, then fill main-sequence holes from the tracks
R=$1; O=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso
F=$O/$R/isos/isochrones_example.txt
[ -s $O/$R/isos/isochrones_pyeep.txt ] && [ -s $F ] && { echo "$R done"; exit 0; }
/dartfs-hpc/scratch/f007058/scmd_tb/build_tuc_12gyr.sh $R > /dev/null 2>&1 || { echo "$R BUILD FAILED"; exit 1; }
mv $F $O/$R/isos/isochrones_pyeep.txt
python /dartfs-hpc/scratch/f007058/scmd_tb/patch_lower_ms.py $O/$R/isos/isochrones_pyeep.txt $O/w$R/mesa $F
rm -rf $O/w$R/raw
