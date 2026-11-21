"""NGC 6791: OG vs new isochrones, one subplot per MC set."""
import sys,os,numpy as np; sys.path.insert(0,'/dartfs-hpc/scratch/f007058/audit_ngc6791')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from trackio import read_mist_iso

A='/dartfs-hpc/scratch/f007058/audit_ngc6791'
OG='/dartfs-hpc/scratch/f007058/subset_isochrones_ngc6791'
runs=[int(x) for x in open(f'{A}/noco100.txt').read().split()]
C_OG,C_NEW='#2a78d6','#eb6834'
plt.rcParams.update({'font.size':8,'axes.grid':True,'grid.alpha':0.2,'grid.linewidth':0.4,
 'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'#fcfcfb',
 'axes.facecolor':'#fcfcfb','savefig.dpi':120})
INK='#0b0b0b'; INK2='#52514e'

def at(path,age):
    cols,a,_=read_mist_iso(path)
    if a is None: return None
    ag=a[:,cols.index('isochrone_age_yr')]; s=a[np.isclose(ag,age,rtol=1e-6)]
    if len(s)==0: return None
    return s[:,cols.index('log_Teff')],s[:,cols.index('log_L')]

for AGE,tag in [(8.0e9,'8gyr'),(12.0e9,'12gyr')]:
    panels=[]
    for r in runs:
        p=f'{A}/pyeep_work/isos/noco{r}/isochrones_noco{r}.txt'
        if not os.path.exists(p): continue
        o=at(f'{OG}/isochrones_{r}.txt',AGE); n=at(p,AGE)
        if o is None: continue
        v={}
        for line in open(f'/dartfs-hpc/rc/home/8/f007058/mcdsep/MC_var_noco/varfeh+031.{r}'):
            q=line.split()
            if len(q)>=2: v[q[1]]=float(q[0])
        panels.append((r,o,n,v))
    N=len(panels); ncol=7; nrow=int(np.ceil(N/ncol))
    fig,axes=plt.subplots(nrow,ncol,figsize=(ncol*2.5,nrow*2.5),sharex=True,sharey=True)
    axes=np.atleast_2d(axes)
    for i,(r,o,n,v) in enumerate(panels):
        ax=axes[i//ncol,i%ncol]
        ax.plot(o[0],o[1],lw=1.3,color=C_OG)
        if n is not None: ax.plot(n[0],n[1],lw=1.3,ls='--',color=C_NEW)
        miss='' if n is not None else '  [no new]'
        ax.set_title('%d  $\\alpha_{MLT}$=%.2f  [Fe/H]=%+.2f%s'%(r,v['CMIXLA'],v['FeH'],miss),
                     fontsize=7,color=(INK if n is not None else '#b00'),pad=2)
        ax.tick_params(labelsize=6,colors=INK2)
    for j in range(N,nrow*ncol): axes[j//ncol,j%ncol].axis('off')
    axes[0,0].invert_xaxis()
    fig.supxlabel('log $T_{\\rm eff}$',color=INK); fig.supylabel('log $L/L_\\odot$',color=INK)
    fig.legend(handles=[Line2D([],[],color=C_OG,lw=2,label='OG isochrone'),
                        Line2D([],[],color=C_NEW,lw=2,ls='--',label='New (overshoot off)')],
               loc='upper right',frameon=False,fontsize=11,ncol=2,bbox_to_anchor=(0.995,0.999))
    fig.suptitle('NGC 6791 — OG vs new isochrones at %.0f Gyr  (%d MC sets)'%(AGE/1e9,N),
                 fontsize=14,color=INK,x=0.005,ha='left',y=0.999)
    fig.tight_layout(rect=[0.012,0.012,1,0.975])
    out=f'{A}/fig10_6791_grid99_{tag}.png'
    fig.savefig(out,bbox_inches='tight'); plt.close(fig)
    print('wrote %s  (%d panels)'%(os.path.basename(out),N))
