import sys,numpy as np; sys.path.insert(0,'/dartfs-hpc/scratch/f007058/audit_ngc6791')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from trackio import read_track
A='/dartfs-hpc/scratch/f007058/audit_ngc6791'
plt.rcParams.update({'font.size':11,'axes.grid':True,'grid.alpha':0.25,'grid.linewidth':0.6,
 'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'#fcfcfb',
 'axes.facecolor':'#fcfcfb','savefig.dpi':150})
INK='#0b0b0b'; INK2='#52514e'
GOOD='#2a78d6'; SICK='#eb6834'; EXT='#1baf7a'

M='1200'
runs=[('var 20000  $\\alpha_C$=0.027  (healthy)',f'{A}/extracted/20000/m{M}.feh031.20000.gz',GOOD,'-'),
      ('var 20006  $\\alpha_C$=0.174  (nmod 9999)',f'{A}/extracted/20006/m{M}.feh031.20006.gz',SICK,'-'),
      ('var 20006  $\\alpha_C$=0.174  (nmod 40000)',f'{A}/nmodtest_x/m{M}.feh031.20006.gz',EXT,'--')]
T=[(lab,read_track(p),c,ls) for lab,p,c,ls in runs]

fig,ax=plt.subplots(2,2,figsize=(13.5,10))

for lab,t,c,ls in T:
    if t is None: continue
    age=t['star_age']/1e9
    ax[0,0].plot(t['log_Teff'],t['log_L'],lw=1.6,color=c,ls=ls,label=lab)
    ax[0,0].plot(t['log_Teff'][-1],t['log_L'][-1],'o',color=c,ms=9,mec='w',mew=1.2)
    ax[0,1].plot(age,t['log_L'],lw=1.6,color=c,ls=ls)
    ax[0,1].plot(age[-1],t['log_L'][-1],'o',color=c,ms=9,mec='w',mew=1.2)
    ax[1,0].plot(age,t['center_h1'],lw=1.6,color=c,ls=ls)
    ax[1,0].plot(age[-1],t['center_h1'][-1],'o',color=c,ms=9,mec='w',mew=1.2)
    dt=np.diff(t['star_age'])
    ax[1,1].plot(np.arange(len(dt)),dt,lw=0.9,color=c,ls=ls,alpha=0.85)

ax[0,0].invert_xaxis(); ax[0,0].set_xlabel('log $T_{\\rm eff}$'); ax[0,0].set_ylabel('log $L/L_\\odot$')
ax[0,0].set_title('A. Sick track leaves the MS but stalls on the subgiant\nbranch, far below the RGB (dot = where it stops)',loc='left',color=INK)
ax[0,0].legend(frameon=False,fontsize=9,loc='upper left')
ax[0,0].axhline(2.6,color=INK2,lw=0.8,ls=':'); ax[0,0].text(0.02,0.90,'RGB tip $\\approx$ 2.6',
    transform=ax[0,0].transAxes,fontsize=9,color=INK2)

ax[0,1].set_xlabel('age (Gyr)'); ax[0,1].set_ylabel('log $L/L_\\odot$')
ax[0,1].set_title('B. Both finish core H burning; only the healthy\none goes on to climb the giant branch',loc='left',color=INK)
ax[0,1].axhline(2.6,color=INK2,lw=0.8,ls=':')

ax[1,0].set_xlabel('age (Gyr)'); ax[1,0].set_ylabel('central hydrogen $X_c$')
ax[1,0].set_title('C. Overshoot extends the MS: TAMS moves\n3.85 $\\rightarrow$ 5.68 Gyr (real, expected physics)',loc='left',color=INK)
ax[1,0].axhline(0,color=INK2,lw=0.8)

ax[1,1].set_yscale('log'); ax[1,1].set_xlabel('model number'); ax[1,1].set_ylabel('timestep $\\Delta t$ (yr)')
ax[1,1].set_title('D. THE CAUSE: at TAMS $\\Delta t$ drops 800x and then\nLOCKS at 1.396e4 yr for every later model',loc='left',color=INK)
ax[1,1].axvline(9999,color=INK2,lw=1.0,ls=':')
ax[1,1].text(10400,2e3,'nmod=9999\ncut-off',fontsize=9,color=INK2)
ax[1,1].axhline(1.396e4,color=SICK,lw=0.8,ls=':')
ax[1,1].text(20000,1.7e4,'locked at 1.396e4 yr',fontsize=9,color=SICK)

for a in ax.ravel(): a.tick_params(colors=INK2)
fig.suptitle('Why high-$\\alpha_C$ tracks truncate — 1.20 $M_\\odot$, NGC 6791',
             fontsize=14,color=INK,x=0.005,ha='left')
fig.tight_layout(rect=[0,0,1,0.95])
fig.savefig(f'{A}/fig5_truncation_mechanism.png',bbox_inches='tight')
print('wrote fig5_truncation_mechanism.png')
for lab,t,c,ls in T:
    if t is None: continue
    dt=np.diff(t['star_age'])
    print('  %-42s models=%6d  end=%5.2f Gyr  maxlogL=%.3f  Xc_end=%.4f  median dt=%.2e'%(
        lab.replace('$\\alpha_C$','aC'),len(t['star_age']),t['star_age'][-1]/1e9,
        t['log_L'].max(),t['center_h1'][-1],np.median(dt)))
