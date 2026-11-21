#!/usr/bin/env python3
"""Truncate each emulated isochrone at ITS OWN RGB tip (max log L per age).

Capping at the real isochrone's max mass does not work: the emulator's
mass->luminosity mapping is offset (the real RGB spans 1.00-1.05 Msun at
8.5 Gyr, the emulator puts it at 1.10-1.20), so a real-mass cap removes the
emulated RGB entirely. Cutting at the emulator's own tip keeps a
self-consistent isochrone -- MS through tip -- and drops only the post-tip
extrapolation into masses that no longer exist.
"""
import os, sys
import numpy as np, pandas as pd

C = ["EEP","isochrone_age_yr","initial_mass","star_mass","star_mdot","he_core_mass",
     "c_core_mass","log_L","log_LH","log_LHe","log_Teff","log_R","log_g","surface_h1",
     "surface_he3","surface_he4","surface_c12","surface_o16","log_center_T",
     "log_center_Rho","center_gamma","center_h1","center_he4","center_c12"]
EMUL = '/dartfs-hpc/scratch/f007058/Isochrones_emul4k'
OUT  = '/dartfs-hpc/scratch/f007058/Isochrones_emul4k_tip'

def main(runs):
    for run in runs:
        src = f'{EMUL}/{run}/isos/isochrones_example.txt'
        if not os.path.exists(src):
            print(f"  {run}: missing"); continue
        dd = f'{OUT}/{run}/isos'; os.makedirs(dd, exist_ok=True)
        dst = f'{dd}/isochrones_example.txt'
        if os.path.exists(dst) and os.path.getsize(dst) > 1e6:
            print(f"  {run}: exists"); continue
        d = pd.read_csv(src, sep=r'\s+', comment='#', skiprows=10); d.columns = C
        keep = []
        for age, s in d.groupby('isochrone_age_yr', sort=False):
            s = s.sort_values('initial_mass')
            tip = int(np.argmax(s.log_L.values))      # the isochrone's own tip
            keep.append(s.iloc[:tip + 1])
        out = pd.concat(keep).reset_index(drop=True)
        hdr = open(src).readlines()[:11]
        npts = len(out) // max(out.isochrone_age_yr.nunique(), 1)
        hdr[8] = f"# number of EEPs, cols = {npts:5d}   24\n"
        with open(dst, 'w') as fh:
            fh.writelines(hdr); np.savetxt(fh, out.values, fmt='%28.16E')
        print(f"  {run}: kept {len(out):,}/{len(d):,} ({len(out)/len(d)*100:.0f}%)  "
              f"{os.path.getsize(dst)/1e6:.0f} MB")

if __name__ == '__main__':
    main(sys.argv[1:])
