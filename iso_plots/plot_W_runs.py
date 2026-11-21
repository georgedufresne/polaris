#!/usr/bin/env python3
"""W (official isochrone_W.py) vs OG / New / old emulator, runs 20000-20002."""
import sys, numpy as np, pandas as pd
sys.path.insert(0,'/dartfs-hpc/rc/home/8/f007058/CONF1DENCE-W')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from isochrone_W import w_isochrone
from predict_W import WEmulator, AgeHead

C=["EEP","isochrone_age_yr","initial_mass","star_mass","star_mdot","he_core_mass",
   "c_core_mass","log_L","log_LH","log_LHe","log_Teff","log_R","log_g","surface_h1",
   "surface_he3","surface_he4","surface_c12","surface_o16","log_center_T",
   "log_center_Rho","center_gamma","center_h1","center_he4","center_c12"]
REAL={'OG':('/dartfs-hpc/rc/lab/C/ChaboyerB/Isochrones_Only_July_10/{r}/isos/isochrones_example.txt','#1f77b4','-'),
      'New':('/dartfs-hpc/scratch/f007058/Isochrones_noco/{r}/isos/isochrones_example.txt','#2ca02c','-'),
      'Emul (old)':('/dartfs-hpc/scratch/f007058/Isochrones_emul4k_tip/{r}/isos/isochrones_example.txt','#ff7f0e','--')}
def ld(p):
    d=pd.read_csv(p,sep=r'\s+',comment='#',skiprows=10); d.columns=C; return d
def sl(d,a):
    x=np.unique(d.isochrone_age_yr); t=x[np.argmin(abs(x-a*1e9))]
    return d[d.isochrone_age_yr==t].sort_values('EEP')

W=WEmulator(); ah=AgeHead(); lh=W.fmt['label_header']
def phys_for(run):
    v=[float(l.split()[0]) for l in open(f'/dartfs-hpc/rc/home/8/f007058/mcdsep/MC_var_noco/varfeh+031.{run}') if l.strip()]
    p=list(W.fmt['label_default'])
    for nm,x in zip(lh[:22],[*v[0:16],1.0,*v[16:21]]): p[lh.index(nm)]=x
    return p

RUNS=['20000','20001','20002']; AGES=[4.0,8.5,12.0]
for tag,zoom in (('full',False),('window',True)):
    fig,axes=plt.subplots(len(RUNS),len(AGES),figsize=(19,16))
    for i,run in enumerate(RUNS):
        ph=phys_for(run)
        for j,age in enumerate(AGES):
            ax=axes[i,j]
            for lab,(tpl,col,ls) in REAL.items():
                try: s=sl(ld(tpl.format(r=run)),age)
                except Exception: continue
                ax.plot(s.log_Teff,s.log_L,ls,color=col,lw=1.5,label=lab)
            iso=w_isochrone(W,ah,ph,age)
            if len(iso['eep']):
                ax.plot(iso['log_Teff'],iso['log_L'],'-',color='#d62728',lw=1.9,label='W')
            if not zoom: ax.invert_xaxis()
            else:
                ax.set_xlim(3.83,3.65); ax.set_ylim(-0.15,1.05)
                ax.axhspan(0.073,0.761,color='gold',alpha=.12,zorder=0)
            ax.grid(alpha=.25); ax.set_title(f'run {run} — {age:.1f} Gyr',fontsize=11)
            ax.set_xlabel('log $T_{eff}$',fontsize=9); ax.set_ylabel('log $L/L_\\odot$',fontsize=9)
            if i==0 and j==0: ax.legend(fontsize=9)
    t='full HRD' if not zoom else 'sCMD fit window (shaded)'
    fig.suptitle(f'NGC 6791 runs 20000-20002 — W (official isochrone_W.py) vs OG / New / old emulator — {t}',fontsize=15)
    fig.tight_layout(rect=[0,0,1,0.97])
    o=f'/dartfs-hpc/scratch/f007058/iso_plots/W_runs20000_20002_{tag}.png'
    fig.savefig(o,dpi=108); plt.close(fig); print('saved',o)

print(f"\n{'run':>7}{'age':>6}{'truth':>8}{'OG':>9}{'New':>9}{'old emul':>10}{'W':>9}{'W err%':>9}")
TRUTH={4.0:None,8.5:None,12.0:None}
for run in RUNS:
    ph=phys_for(run)
    for age in AGES:
        row=f"{run:>7}{age:>6.1f}"
        vals={}
        for lab,(tpl,_,_) in REAL.items():
            try: vals[lab]=sl(ld(tpl.format(r=run)),age).initial_mass.max()
            except Exception: vals[lab]=np.nan
        iso=w_isochrone(W,ah,ph,age); wm=iso['mass'].max() if len(iso['mass']) else np.nan
        truth=vals['New']   # New reproduces the tracks to 3 decimals
        row+=f"{truth:>8.3f}{vals['OG']:>9.3f}{vals['New']:>9.3f}{vals['Emul (old)']:>10.3f}{wm:>9.3f}{(wm/truth-1)*100:>9.1f}"
        print(row)
