#!/usr/bin/env python3
"""Emulated vs real NGC 6791 isochrones."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
C=["EEP","isochrone_age_yr","initial_mass","star_mass","star_mdot","he_core_mass",
   "c_core_mass","log_L","log_LH","log_LHe","log_Teff","log_R","log_g","surface_h1",
   "surface_he3","surface_he4","surface_c12","surface_o16","log_center_T",
   "log_center_Rho","center_gamma","center_h1","center_he4","center_c12"]
def load(p):
    d=pd.read_csv(p,sep=r'\s+',comment='#',skiprows=10); d.columns=C; return d
SRC={
 'OG':      ('/dartfs-hpc/rc/lab/C/ChaboyerB/Isochrones_Only_July_10/{r}/isos/isochrones_example.txt','#1f77b4'),
 'New':     ('/dartfs-hpc/scratch/f007058/Isochrones_noco/{r}/isos/isochrones_example.txt','#2ca02c'),
 'Emul 4k': ('/dartfs-hpc/scratch/f007058/Isochrones_emul4k/{r}/isos/isochrones_example.txt','#d62728'),
}
RUN='20000'; AGES=[5.0,8.5,12.0]
fig,axes=plt.subplots(2,3,figsize=(19,11))
for j,age in enumerate(AGES):
    top,bot=axes[0,j],axes[1,j]
    for lab,(tpl,col) in SRC.items():
        try: d=load(tpl.format(r=RUN))
        except Exception: continue
        a=np.unique(d.isochrone_age_yr); t=a[np.argmin(abs(a-age*1e9))]
        s=d[d.isochrone_age_yr==t].sort_values('initial_mass')
        top.plot(s.log_Teff,s.log_L,'-',color=col,lw=1.5,
                 label=f"{lab} ({len(s)} pts, m<{s.initial_mass.max():.2f})")
        bot.plot(s.log_Teff,s.log_L,'-',color=col,lw=1.8,label=lab)
    top.set_title(f'{age:.1f} Gyr — full'); top.invert_xaxis()
    bot.set_title(f'{age:.1f} Gyr — sCMD fit window')
    bot.set_xlim(3.82,3.66); bot.set_ylim(-0.1,1.0)
    bot.axhspan(0.073,0.761,color='gold',alpha=.13,zorder=0)
    for ax in (top,bot):
        ax.grid(alpha=.25); ax.set_xlabel('log $T_{eff}$'); ax.set_ylabel('log $L/L_\\odot$')
        ax.legend(fontsize=8)
fig.suptitle(f'NGC 6791 run {RUN}: emulated (central, 4000 pts) vs real isochrones',fontsize=15)
fig.tight_layout(rect=[0,0,1,0.96])
o='/dartfs-hpc/scratch/f007058/emul_iso/iso_emul_vs_real_20000.png'
fig.savefig(o,dpi=110); print('saved',o)
# mass extent table
print(f"\n{'age':>6}  " + "  ".join(f"{k:>22s}" for k in SRC))
for age in [4.0,6.0,8.5,10.0,12.0]:
    row=f"{age:>6.1f}  "
    for lab,(tpl,_) in SRC.items():
        try:
            d=load(tpl.format(r=RUN)); a=np.unique(d.isochrone_age_yr)
            s=d[d.isochrone_age_yr==a[np.argmin(abs(a-age*1e9))]]
            row+=f"{('m≤%.3f n=%d'%(s.initial_mass.max(),len(s))):>22s}  "
        except Exception: row+=f"{'-':>22s}  "
    print(row)
