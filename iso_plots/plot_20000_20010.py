#!/usr/bin/env python3
"""OG vs New vs Emulated isochrones, runs 20000-20010."""
import numpy as np, pandas as pd, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
C=["EEP","isochrone_age_yr","initial_mass","star_mass","star_mdot","he_core_mass",
   "c_core_mass","log_L","log_LH","log_LHe","log_Teff","log_R","log_g","surface_h1",
   "surface_he3","surface_he4","surface_c12","surface_o16","log_center_T",
   "log_center_Rho","center_gamma","center_h1","center_he4","center_c12"]
SRC={'OG':('/dartfs-hpc/rc/lab/C/ChaboyerB/Isochrones_Only_July_10/{r}/isos/isochrones_example.txt','#1f77b4'),
     'New':('/dartfs-hpc/scratch/f007058/Isochrones_noco/{r}/isos/isochrones_example.txt','#2ca02c'),
     'Emulated':('/dartfs-hpc/scratch/f007058/Isochrones_emul4k_tip/{r}/isos/isochrones_example.txt','#d62728')}
OUT='/dartfs-hpc/scratch/f007058/iso_plots'
RUNS=[str(r) for r in range(20000,20011)]
WLO,WHI=0.073,0.761
def load(p):
    d=pd.read_csv(p,sep=r'\s+',comment='#',skiprows=10); d.columns=C; return d
def sl(d,age):
    a=np.unique(d.isochrone_age_yr); return d[d.isochrone_age_yr==a[np.argmin(abs(a-age*1e9))]].sort_values('initial_mass')

for AGE in (8.5,):
    for tag,fn in [('full',f'runs20000_20010_{AGE}Gyr_full.png'),
                   ('window',f'runs20000_20010_{AGE}Gyr_window.png')]:
        fig,axes=plt.subplots(3,4,figsize=(21,14))
        for ax,r in zip(axes.ravel(),RUNS):
            for lab,(tpl,col) in SRC.items():
                try: s=sl(load(tpl.format(r=r)),AGE)
                except Exception: continue
                ax.plot(s.log_Teff,s.log_L,'-',color=col,lw=1.5,label=lab)
            ax.set_title(f'run {r}',fontsize=11); ax.grid(alpha=.25)
            if tag=='full': ax.invert_xaxis()
            else:
                ax.set_xlim(3.82,3.66); ax.set_ylim(-0.1,1.0)
                ax.axhspan(WLO,WHI,color='gold',alpha=.12,zorder=0)
            ax.set_xlabel('log $T_{eff}$',fontsize=9); ax.set_ylabel('log $L/L_\\odot$',fontsize=9)
            ax.tick_params(labelsize=8)
        for ax in axes.ravel()[len(RUNS):]: ax.axis('off')
        axes.ravel()[0].legend(fontsize=10)
        t='full HRD' if tag=='full' else 'sCMD fit window (shaded)'
        fig.suptitle(f'NGC 6791 isochrones, runs 20000-20010 at {AGE} Gyr — OG vs New vs Emulated — {t}',fontsize=16)
        fig.tight_layout(rect=[0,0,1,0.965]); fig.savefig(f'{OUT}/{fn}',dpi=105); plt.close(fig)
        print('saved',fn)

# per-run numeric offsets in the fit window
grid=np.linspace(WLO,WHI,80); rows=[]
for r in RUNS:
    try:
        cur={}
        for lab,(tpl,_) in SRC.items():
            s=sl(load(tpl.format(r=r)),8.5).drop_duplicates('log_L').sort_values('log_L')
            cur[lab]=np.interp(grid,s.log_L,s.log_Teff,left=np.nan,right=np.nan)
        rows.append({'run':r,
                     'New-OG':np.nanmean(cur['New']-cur['OG']),
                     'Emul-OG':np.nanmean(cur['Emulated']-cur['OG']),
                     'Emul-New':np.nanmean(cur['Emulated']-cur['New'])})
    except Exception as e: print(f'  {r}: {e}')
df=pd.DataFrame(rows)
print(f"\nmean Δlog Teff in the fit window at 8.5 Gyr (ΔT at 5500 K in brackets)")
print(f"{'run':>7}{'New-OG':>20}{'Emul-OG':>20}{'Emul-New':>20}")
for _,x in df.iterrows():
    print(f"{x.run:>7}" + "".join(f"{f'{x[k]:+.5f} ({5500*(10**x[k]-1):+.0f} K)':>20}" for k in ['New-OG','Emul-OG','Emul-New']))
print(f"{'MEAN':>7}" + "".join(f"{f'{df[k].mean():+.5f} ({5500*(10**df[k].mean()-1):+.0f} K)':>20}" for k in ['New-OG','Emul-OG','Emul-New']))
df.to_csv(f'{OUT}/offsets_20000_20010.csv',index=False)
