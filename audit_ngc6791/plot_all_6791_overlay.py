"""All NGC 6791 isochrones on one plot: every OG vs every new, 8 and 12 Gyr."""
import sys,os,numpy as np; sys.path.insert(0,'/dartfs-hpc/scratch/f007058/audit_ngc6791')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from trackio import read_mist_iso

A='/dartfs-hpc/scratch/f007058/audit_ngc6791'
OG='/dartfs-hpc/scratch/f007058/subset_isochrones_ngc6791'
runs=[int(x) for x in open(f'{A}/noco50.txt').read().split()]
C_OG,C_NEW='#2a78d6','#eb6834'
plt.rcParams.update({'font.size':12,'axes.grid':True,'grid.alpha':0.25,'grid.linewidth':0.6,
 'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'#fcfcfb',
 'axes.facecolor':'#fcfcfb','savefig.dpi':150})
INK='#0b0b0b'; INK2='#52514e'

def at(path,age):
    cols,a,_=read_mist_iso(path)
    if a is None: return None
    ag=a[:,cols.index('isochrone_age_yr')]
    s=a[np.isclose(ag,age,rtol=1e-6)]
    if len(s)==0: return None
    return s[:,cols.index('log_Teff')],s[:,cols.index('log_L')]

fig,ax=plt.subplots(1,2,figsize=(15,6.8))
nn=0
for j,age in enumerate([8.0e9,12.0e9]):
    k=0
    for r in runs:
        p=f'{A}/pyeep_work/isos/noco{r}/isochrones_noco{r}.txt'
        if not os.path.exists(p): continue
        o=at(f'{OG}/isochrones_{r}.txt',age); n=at(p,age)
        if o is None or n is None: continue
        ax[j].plot(o[0],o[1],lw=0.8,color=C_OG,alpha=0.45)
        ax[j].plot(n[0],n[1],lw=0.8,color=C_NEW,alpha=0.45,ls='--')
        k+=1
    nn=max(nn,k)
    ax[j].invert_xaxis()
    ax[j].set_xlabel('log $T_{\\rm eff}$'); ax[j].set_ylabel('log $L/L_\\odot$')
    ax[j].set_title('%.0f Gyr  —  %d MC sets overlaid'%(age/1e9,k),loc='left',color=INK)
    ax[j].tick_params(colors=INK2)
ax[0].legend(handles=[Line2D([],[],color=C_OG,lw=2,label='OG isochrones'),
                      Line2D([],[],color=C_NEW,lw=2,ls='--',label='New (overshoot off)')],
             frameon=False,loc='upper left',fontsize=11)
fig.suptitle('NGC 6791 — all isochrones overlaid, OG vs new',fontsize=15,color=INK,x=0.005,ha='left')
fig.tight_layout(rect=[0,0,1,0.94])
out=f'{A}/fig7_all_6791_overlay.png'
fig.savefig(out,bbox_inches='tight'); print('wrote fig7_all_6791_overlay.png  (%d sets)'%nn)
