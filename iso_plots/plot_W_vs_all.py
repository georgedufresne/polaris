#!/usr/bin/env python3
"""W vs the old emulator vs OG vs New, NGC 6791 run 20000."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
C=["EEP","isochrone_age_yr","initial_mass","star_mass","star_mdot","he_core_mass",
   "c_core_mass","log_L","log_LH","log_LHe","log_Teff","log_R","log_g","surface_h1",
   "surface_he3","surface_he4","surface_c12","surface_o16","log_center_T",
   "log_center_Rho","center_gamma","center_h1","center_he4","center_c12"]
S={'OG':('/dartfs-hpc/rc/lab/C/ChaboyerB/Isochrones_Only_July_10/{r}/isos/isochrones_example.txt','#1f77b4','-'),
   'New':('/dartfs-hpc/scratch/f007058/Isochrones_noco/{r}/isos/isochrones_example.txt','#2ca02c','-'),
   'Emul (old, age-cond)':('/dartfs-hpc/scratch/f007058/Isochrones_emul4k_tip/{r}/isos/isochrones_example.txt','#ff7f0e','--'),
   'W (EEP-cond)':('/dartfs-hpc/scratch/f007058/Isochrones_W/{r}/isos/isochrones_example.txt','#d62728','-')}
def ld(p):
    d=pd.read_csv(p,sep=r'\s+',comment='#',skiprows=10); d.columns=C; return d
def sl(d,a):
    x=np.unique(d.isochrone_age_yr); t=x[np.argmin(abs(x-a*1e9))]
    return d[d.isochrone_age_yr==t].sort_values('EEP')
RUN='20000'; WLO,WHI=0.073,0.761
fig,axes=plt.subplots(2,3,figsize=(20,11))
for j,age in enumerate([5.0,8.5,12.0]):
    for i,(ax,zoom) in enumerate([(axes[0,j],False),(axes[1,j],True)]):
        for lab,(tpl,col,ls) in S.items():
            try: s=sl(ld(tpl.format(r=RUN)),age)
            except Exception: continue
            ax.plot(s.log_Teff,s.log_L,ls,color=col,lw=1.6,
                    label=f"{lab} (m≤{s.initial_mass.max():.3f})" if not zoom else lab)
        if not zoom:
            ax.invert_xaxis(); ax.set_title(f'{age:.1f} Gyr — full')
        else:
            ax.set_xlim(3.82,3.66); ax.set_ylim(-0.1,1.0)
            ax.axhspan(WLO,WHI,color='gold',alpha=.12,zorder=0)
            ax.set_title(f'{age:.1f} Gyr — sCMD fit window')
        ax.grid(alpha=.25); ax.set_xlabel('log $T_{eff}$'); ax.set_ylabel('log $L/L_\\odot$')
        ax.legend(fontsize=8)
fig.suptitle('NGC 6791 run 20000 — W (EEP-conditioned) vs the old age-conditioned emulator vs OG / New',fontsize=15)
fig.tight_layout(rect=[0,0,1,0.955])
o='/dartfs-hpc/scratch/f007058/iso_plots/W_vs_old_vs_OG_20000.png'
fig.savefig(o,dpi=110); print('saved',o)
# numbers
print(f"\n{'age':>6}  " + "".join(f"{k:>26s}" for k in S))
for age in [4.0,6.0,8.5,10.0,12.0]:
    row=f"{age:>6.1f}  "
    for lab,(tpl,_,_) in S.items():
        try:
            s=sl(ld(tpl.format(r=RUN)),age)
            row+=f"{f'm<={s.initial_mass.max():.3f} n={len(s)}':>26s}"
        except Exception: row+=f"{'-':>26s}"
    print(row)
print("\ntrack truth max alive mass: 4.0 Gyr 1.296 | 8.5 Gyr 1.052 | 12.0 Gyr 0.957")
