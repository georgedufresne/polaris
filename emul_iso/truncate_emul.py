#!/usr/bin/env python3
"""Cap each emulated isochrone's mass range, per age, at the real isochrone's.

Why: scmd_grid_...py samples the IMF over the isochrone's OWN mass range
(m_max = np.max(initial_mass), line ~302). The emulator is asked for a fixed
0.25-1.60 Msun grid at every age, but a real NGC 6791 isochrone stops at the
RGB tip -- 1.05 Msun at 8.5 Gyr, 0.96 at 12 Gyr. Everything above that is
extrapolation into masses that are already dead, and it has plausible log L so
the value filter cannot catch it. Left in, the sCMD populates the CMD with
stars that should not exist and the 2D KS lands at ~0.60 instead of ~0.10.

Truncating at the real per-run, per-age maximum mass makes the emulated and
real sCMDs sample the same mass range, so the comparison measures the emulator
rather than this artefact. Note this does use the real grid to set the extent:
the emulator cannot self-truncate, which is itself worth reporting.
"""
import glob, json, os, sys
import numpy as np, pandas as pd

C = ["EEP","isochrone_age_yr","initial_mass","star_mass","star_mdot","he_core_mass",
     "c_core_mass","log_L","log_LH","log_LHe","log_Teff","log_R","log_g","surface_h1",
     "surface_he3","surface_he4","surface_c12","surface_o16","log_center_T",
     "log_center_Rho","center_gamma","center_h1","center_he4","center_c12"]
EMUL = '/dartfs-hpc/scratch/f007058/Isochrones_emul4k'
REAL = '/dartfs-hpc/scratch/f007058/Isochrones_noco'
OUT  = '/dartfs-hpc/scratch/f007058/Isochrones_emul4k_cut'
MED  = json.load(open('/dartfs-hpc/scratch/f007058/emul_iso/massmax_vs_age.json'))
MED  = {float(k): v for k, v in MED.items()}

def read(p):
    d = pd.read_csv(p, sep=r'\s+', comment='#', skiprows=10)
    d.columns = C
    return d

def cap_for(run):
    """per-age max mass from this run's own real isochrone; median fallback."""
    p = f'{REAL}/{run}/isos/isochrones_example.txt'
    if os.path.exists(p):
        try:
            r = read(p)
            return r.groupby(r.isochrone_age_yr.round(-7))['initial_mass'].max().to_dict(), 'own'
        except Exception:
            pass
    return MED, 'median'

def main(runs):
    for run in runs:
        src = f'{EMUL}/{run}/isos/isochrones_example.txt'
        if not os.path.exists(src):
            print(f"  {run}: missing"); continue
        dst_d = f'{OUT}/{run}/isos'; os.makedirs(dst_d, exist_ok=True)
        dst = f'{dst_d}/isochrones_example.txt'
        if os.path.exists(dst) and os.path.getsize(dst) > 1e6:
            print(f"  {run}: exists"); continue
        caps, how = cap_for(run)
        d = read(src)
        keyed = d.isochrone_age_yr.round(-7)
        cap_series = keyed.map(lambda a: caps.get(a, MED.get(a, 1.6)))
        keep = d.initial_mass.values <= cap_series.values
        out = d[keep].copy()
        hdr = open(src).readlines()[:11]
        npts = len(out) // max(out.isochrone_age_yr.nunique(), 1)
        hdr[8] = f"# number of EEPs, cols = {npts:5d}   24\n"
        with open(dst, 'w') as fh:
            fh.writelines(hdr)
            np.savetxt(fh, out.values, fmt='%28.16E')
        print(f"  {run}: kept {len(out):,}/{len(d):,} ({len(out)/len(d)*100:.0f}%) "
              f"cap={how}  m_max {out.initial_mass.max():.3f}  {os.path.getsize(dst)/1e6:.0f} MB")

if __name__ == '__main__':
    main(sys.argv[1:])
