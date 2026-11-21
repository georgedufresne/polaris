#! /bin/bash 
#
# script which automatically runs the Dartmouth stellar evolution code  
# for a Monte Carlo simulation
# requires program newpoly 
# requires subsrunmc.sh 
# requires gs98getz to change [Fe/H] & [alpha/Fe]
#

#export LD_LIBRARY_PATH=/dartfs-hpc/admin/users/chaboyer/free_eos-2.2.1/intelbuild:$LD_LIBRARY_PATH
#export LD_LIBRARY_PATH=/dartfs-hpc/rc/lab/C/ChaboyerB/Catherine/dsep3/free_eos-2.2.1/intelbuild:$LD_LIBRARY_PATH
#LD_LIBRARY_PATH=/dartfs-hpc/rc/home/v/f0056bv/intelbuild:$LD_LIBRARY_PATH
#LD_LIBRARY_PATH=/afs/northstar/users/c/chaboyer/public/intelbuild:$LD_LIBRARY_PATH
#LD_LIBRARY_PATH=/home1/08819/mying/dsepX-master/DSEPmisc/libs/free_eos-2.2.1/build/:$LD_LIBRARY_PATH
# base directory; most other directories are below this
#local=/dartfs-hpc/rc/home/v/f0056bv/dsep3
local=/dartfs-hpc/scratch/f007058/M15/dsep4/env/DSEPmisc

poly=$local/poly
zams=$local/zams
prems=$local/prems
opac=$local/opac
atm=$local/surfBC
opal=$opac/OPAL
phx=$opac/phoenix
bcphx=$atm/phoenix/GS98
bckur=$atm/atmk

GC_name=47Tuc

#mc directory
# NML can point at a retry namelist set (e.g. damped first corrections)
nml=${NML:-/dartfs-hpc/scratch/f007058/M15/dsep4/env/nml}
#nml=/Users/george/DSEPmisc/nml
#/Users/george/DSEPmisc/nml
#nml=/Users/george/Documents/MCDSEP/dsep3/nml
#nml=$local/nml
prog=/dartfs-hpc/scratch/f007058/M15/dsep4/env/dsep4/bin
#vardir=/dartfs-hpc/rc/lab/C/ChaboyerB/Catherine/var
#vardir=/OneDrive/Desktop/M55
vardir=/dartfs-hpc/scratch/f007058/scmd_tb/tuc_d4v2_build/var
rundir=/dartfs-hpc/scratch/f007058/M15/dsep4/env/run

#read from command line which run to do
if [ $# -ne 3 ]; then
    echo "Usage: $0 <varnum> <feh> <mass>"
    exit 1
fi

# create a  random working directory

number=$RANDOM
num2=$(date +%s)
num2=$((num2 - 1265900000))

# unique per task: $RANDOM+seconds can collide when hundreds of array tasks start together
out=$(mktemp -d /dartfs-hpc/scratch/f007058/tuc_dsep4_batch/work/d4fix_$1_$3_XXXXXX)
echo "Scratch dir: $out"
# (mktemp already created $out)

# directory which output gets copied to from the temp directory
#storeoutput=/dartfs-hpc/rc/home/8/f007058/mcdsep/$GC_name/tracks
storeoutput=${STOREBASE:-/dartfs-hpc/scratch/f007058/Tracks_tuc_dsep4_mcfix}/$1
#/dartfs-hpc/rc/lab/C/ChaboyerB/NGC6791/tracks
echo "Output dir: $storeoutput"

ID=GS98
Lmass=($(seq 200 40 690))
MLmass=($(seq 700 50 1400))
MHmass=($(seq 1500 100 1900))
Hmass=($(seq 2000 200 3000))
masses=( $3 )   # one mass per task: $3
#Run quick test
#masses=($(seq 1000 50 1001))
nmod=60000
endage=1.8D10

runum=$1
tempfeh=$2
sign=${tempfeh:0:1}
absfeh=${tempfeh:1}
if [[ $sign != '-' ]]; then
  sign="+"
  absfeh=$tempfeh
fi
intfeh=${absfeh//./}
#tmpfeh=${2#-}
#absfeh=`echo "$tmpfeh * 100" | bc`
#intfeh=${absfeh%.*}
#intfeh=77
#varfile=$vardir/varfeh-$intfeh.$runum
#varfile=$varfeh-190.
#varfile=$vardir/varfeh-$intfeh.$runum
varfile=$vardir/varfeh${sign}${intfeh}.$runum
echo "$varfile"

feh=$(cat $varfile | awk 'NR==1 {print $1}')
echo "[Fe/H] = $feh"
afe=$(cat $varfile | awk 'NR==3 {print $1}')
echo "[a/Fe] = $afe"
yy=$(cat $varfile | awk 'NR==2 {print $1}')
echo "Yp = $yy"
mixlen=$(cat $varfile | awk 'NR==4 {print $1}')
echo "mixL = $mixlen"

zz=$($rundir/gs98getz $feh $afe)
echo "Z= $zz"

xx=$(awk "BEGIN {print 1.0-$zz-$yy}")
echo "X= $xx"

#check if file exist. If so, stop the code from running
mkdir -p $storeoutput
if (( $3 < 1000 )); then _mp=0$3; else _mp=$3; fi
# a failed run still leaves a 39-byte gzip of an empty file, so test the content
if [ -n "$(zcat $storeoutput/m$_mp.feh$intfeh.$runum.gz 2>/dev/null | head -c 1)" ]; then
  echo "mass $3 already done for $runum." && exit 0
fi

#check if isochrone is already exist
#iso_dir=/dartfs-hpc/rc/lab/C/ChaboyerB/Catherine/47Tuc/outiso/feh077.$runum
#if [ -f $iso_dir ]; then
#  echo "iso exists." && exit 0
#else
#  echo "iso does not exist and DSEP is about to start"
#fi



#input various subroutines

. $rundir/subsrunmc.sh


cp $rundir/newpoly $out/.
cd $out

#setup a bunch of directories, filenames and parameters for newpoly
setups

# dsep4 Blouin20 conductive opacity: THREE tables (wd/sd/I) on units 75/76/77,
# unlike dsepX_fixed which had a single table on unit 77.
D4OPAC=/dartfs-hpc/scratch/f007058/M15/dsep4/env/dsep4/opac
ln -sf $D4OPAC/condtab21wd.dat fort.75
ln -sf $D4OPAC/condtab21sd.dat fort.76
ln -sf $D4OPAC/condtab21_I.dat fort.77

#setup mixture specific variables
setupmix

#calculate adjustment to starting Teff & Lum based upon Z of the model
# note that Z is calculated in the subroutine setupmix
#bash can't do real number arithmetic, so have to use awk

# These are the originals
#Ladjust=$(awk "BEGIN {print -0.115*log($zz/$zsolar)/log(10.0)}")
#Tadjust=$(awk "BEGIN {print -0.022*log($zz/$zsolar)/log(10.0)}")

# These are Karabo's and they don't work for high metallicity
#Ladjust=$(awk "BEGIN {print -0.015*log($zz/$zsolar)/log(10.0)}")
#Tadjust=$(awk "BEGIN {print -0.010*log($zz/$zsolar)/log(10.0)}")
#Tadjust=$(awk "BEGIN {print -0.010*log($zz/$zsolar)/log(10.0) + 0.1*($mixlen-1.9)}")

# These are mine and I'm not sure what they are for...
#Ladjust=$(awk "BEGIN {print 0.13141*$mixlen+1.4}")
#Tadjust=$(awk "BEGIN {print 0.04724*$mixlen+3.55}")

#Ladjust=$(awk "BEGIN {print 0.13141*$mixlen+1.35}") #Change Feb-11-2025
#Tadjust=$(awk "BEGIN {print 0.04724*$mixlen+3.50}") #Change Feb-11-2025
#Ladjust=$(awk "BEGIN {print 0.13141*$mixlen+1.36}") #Change Feb-19-2025
#Tadjust=$(awk "BEGIN {print 0.04724*$mixlen+3.51}") #Change Feb-19-2025
#Ladjust=$(awk "BEGIN {print 0.13141*$mixlen+1.37}") #Change March-1-2025
#Tadjust=$(awk "BEGIN {print 0.04724*$mixlen+3.52}") #Change March-1-2025
Ladjust=$(awk "BEGIN {print 0.13141*$mixlen+1.375}") #Change March-2-2025
Tadjust=$(awk "BEGIN {print 0.04724*$mixlen+3.525}") #Change March-2-2025

#Original above

#Ladjust=$(awk "BEGIN {print 0.12*$mixlen+1.5197}")
#Tadjust=$(awk "BEGIN {print 0.042*$mixlen+3.5936}")


#
#Ladjust=$(awk "BEGIN {print 0.22632*$mixlen+0.66413}")
#Tadjust=$(awk "BEGIN {print 0.060983*$mixlen+3.4889}")
#Ladjust=$(awk "BEGIN {print 0.11703*$mixlen+1.0062}")
#Tadjust=$(awk "BEGIN {print 0.073084*$mixlen+3.458}")

ln -s $varfile    fort.85


for mass in "${masses[@]}"; do

echo "mass = $mass"

#calculate Teff, and luminosity  of the starting pre-MS model
    sumass=$(awk "BEGIN {print $mass/1000.0}")
    get_teff_lum
    # START=massdep: mass-dependent pre-MS seed (the formula get_teff_lum computes
    # and then overrides with the mixing-length-only Ladjust/Tadjust). Used as a
    # retry when DSEP4 cannot converge on the standard seed.
    if [ "${START:-}" = "massdep" ]; then
      if (( mass <= 390 )); then
        teff=$(awk "BEGIN {print 0.430 * $sumass + 3.410}")
        lum=$(awk "BEGIN {print  6.20 * $sumass  -1.20}")
      else
        teff=$(awk "BEGIN {print 3.88 + (1.0 - $sumass)*0.05}")
        lum=$(awk "BEGIN {print 2.9 + 0.2 * ($sumass-5.0)}")
      fi
    fi
    echo "Mass= $sumass  lumonisty=$lum  Teff=$teff"

#write out the namelist and  run the program
    wrtpoly
    echo "Finished wrtpoly"
    ./newpoly
    rm poly.nml

#run stellar evolution code

#assign physics namelist based upon mass
    getphysnml

    if (($mass < 1000 )); then
	mprefix=0$mass
    else
	mprefix=$mass
    fi


    fname=m$mprefix.feh$intfeh.$runum
    fnameiso=iso_m$mprefix.feh$intfeh.$runum

    ln -s $out/$fname fort.35
    ln -s $out/$fnameiso fort.37

    echo "about to start DSEP"
    time $prog/dsep
    rc_dsep=$?
    echo "dsep exit status $rc_dsep"
    gzip $out/$fname
    gzip $out/$fnameiso

#    mv $out/$fname.gz $storeoutput/.
 
    rm fort.{11,12,13,35,19,20,22,37}
 
done
 
cd $out
 
# exit 0 = clean run; 3 = dsep crashed but left a partial track (copied);
# 1 = no track at all (nothing copied, so a retry is not skipped)
if [ -z "$(zcat $out/m$_mp.feh$intfeh.$runum.gz 2>/dev/null | head -c 1)" ]; then
    echo "ERROR: no track for $runum mass $3 (dsep exit $rc_dsep)" >&2
    rm -rf $out
    exit 1
fi
cp -f $out/m$_mp.feh$intfeh.$runum.gz     $storeoutput/.
cp -f $out/iso_m$_mp.feh$intfeh.$runum.gz $storeoutput/. 2>/dev/null
rm -rf $out
[ "$rc_dsep" -eq 0 ] || exit 3
exit 0
