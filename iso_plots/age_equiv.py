#!/usr/bin/env python3
"""Which REAL age does an emulated isochrone of age A actually look like?

If the emulator's isochrone labelled 8.5 Gyr resembles the real isochrone at,
say, 8.75 Gyr, then a CMD fit will recover a YOUNGER age for the same data --
which would explain a shift in the posterior with no visible offset in the
isochrone locus at matched age.

Distance is measured in the sCMD fit window only (log L 0.073-0.761), as
mean |delta log Teff| on a common log L ladder.
"""
import numpy as np, pandas as pd
C=["EEP","isochrone_age_yr","initial_mass","star_mass","star_mdot","he_core_mass",
   "c_core_mass","log_L","log_LH","log_LHe","log_Teff","log_R","log_g","surface_h1",
   "surface_he3","surface_he4","surface_c12","surface_o16","log_center_T",
   "log_center_Rho","center_gamma","center_h1","center_he4","center_c12"]
REAL='/dartfs-hpc/scratch/f007058/Isochrones_noco/{r}/isos/isochrones_example.txt'
EMUL='/dartfs-hpc/scratch/f007058/Isochrones_emul4k_tip/{r}/isos/isochrones_example.txt'
WLO,WHI=0.073,0.761; GRID=np.linspace(WLO,WHI,60)

def load(p):
    d=pd.read_csv(p,sep=r'\s+',comment='#',skiprows=10); d.columns=C; return d
def curve(d,age):
    a=np.unique(d.isochrone_age_yr); t=a[np.argmin(abs(a-age*1e9))]
    s=d[d.isochrone_age_yr==t].drop_duplicates('log_L').sort_values('log_L')
    return np.interp(GRID,s.log_L,s.log_Teff,left=np.nan,right=np.nan)

rows=[]
for run in range(20000,20011):
    try:
        R,E=load(REAL.format(r=run)),load(EMUL.format(r=run))
    except Exception: continue
    real_ages=np.unique(R.isochrone_age_yr)/1e9
    for A in [6.0,7.0,8.0,8.5,9.0,10.0,11.0]:
        e=curve(E,A)
        if np.isnan(e).all(): continue
        best,bd=np.nan,np.inf
        for ra in real_ages:
            r=curve(R,ra)
            m=~(np.isnan(e)|np.isnan(r))
            if m.sum()<20: continue
            dd=np.mean(np.abs(e[m]-r[m]))
            if dd<bd: bd,best=dd,ra
        rows.append({'run':run,'emul_age':A,'best_real_age':best,
                     'equiv_shift':best-A,'dist':bd})
df=pd.DataFrame(rows)
df.to_csv('/dartfs-hpc/scratch/f007058/iso_plots/age_equivalence.csv',index=False)
print("Emulated isochrone labelled age A best matches the REAL isochrone at age B")
print("(fit window only; positive shift = emulator looks OLDER than its label)\n")
print(f"{'emul age':>10}{'n':>4}{'median best real age':>23}{'median shift':>15}{'mean |dlogTeff|':>18}")
for A,g in df.groupby('emul_age'):
    print(f"{A:>10.1f}{len(g):>4}{g.best_real_age.median():>23.2f}{g.equiv_shift.median():>+15.2f}{g.dist.median():>18.5f}")
print(f"\noverall median equivalent-age shift: {df.equiv_shift.median():+.3f} Gyr")
print(f"   mean: {df.equiv_shift.mean():+.3f} +/- {df.equiv_shift.std():.3f} Gyr")
