import sys,os,numpy as np; sys.path.insert(0,'/dartfs-hpc/scratch/f007058/audit_ngc6791')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from trackio import read_mist_iso, read_dsep_iso

A='/dartfs-hpc/scratch/f007058/audit_ngc6791'
AGE_MYR=8000.0
runs=[int(x) for x in open(f'{A}/tuc_good50.txt').read().split()]
C_OG,C_NEW='#2a78d6','#eb6834'
plt.rcParams.update({'font.size':8,'axes.grid':True,'grid.alpha':0.2,'grid.linewidth':0.4,
 'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'#fcfcfb',
 'axes.facecolor':'#fcfcfb','savefig.dpi':120})
INK='#0b0b0b'; INK2='#52514e'

panels=[]
for r in runs:
    d=read_dsep_iso(f'/dartfs-hpc/scratch/f007058/47Tuc_iso/feh077.{r}')
    if AGE_MYR not in d: continue
    o=d[AGE_MYR]
    p=f'{A}/pyeep_work/isos/tuc{r}/isochrones_tuc{r}.txt'
    if not os.path.exists(p): continue
    cols,a,_=read_mist_iso(p)
    ag=a[:,cols.index('isochrone_age_yr')]
    s=a[np.isclose(ag,AGE_MYR*1e6,rtol=1e-6)]
    if len(s)==0: continue
    v={}
    for line in open(f'/dartfs-hpc/scratch/f007058/47Tuc_var/varfeh-077.{r}'):
        q=line.split()
        if len(q)>=2: v[q[1]]=float(q[0])
    panels.append((r,o,s,cols,v))
print('panels:',len(panels))

n=len(panels); ncol=7; nrow=int(np.ceil(n/ncol))
fig,axes=plt.subplots(nrow,ncol,figsize=(ncol*2.5,nrow*2.5),sharex=True,sharey=True)
axes=np.atleast_2d(axes)
for i,(r,o,s,cols,v) in enumerate(panels):
    ax=axes[i//ncol,i%ncol]
    ax.plot(o['Log_T'],o['Log_L'],lw=1.3,color=C_OG,label='OG')
    ax.plot(s[:,cols.index('log_Teff')],s[:,cols.index('log_L')],lw=1.3,ls='--',color=C_NEW,label='new')
    ax.set_title('%d   $\\alpha_{MLT}$=%.2f  [Fe/H]=%+.2f'%(r,v['CMIXLA'],v['FeH']),
                 fontsize=7,color=INK,pad=2)
    ax.tick_params(labelsize=6,colors=INK2)
for j in range(n,nrow*ncol): axes[j//ncol,j%ncol].axis('off')
axes[0,0].invert_xaxis()
fig.supxlabel('log $T_{\\rm eff}$',color=INK); fig.supylabel('log $L/L_\\odot$',color=INK)
fig.legend(handles=[Line2D([],[],color=C_OG,lw=2,label='OG isochrone'),
                    Line2D([],[],color=C_NEW,lw=2,ls='--',label='New isochrone (overshoot off)')],
           loc='upper right',frameon=False,fontsize=11,ncol=2,bbox_to_anchor=(0.995,0.999))
fig.suptitle('47 Tuc — OG vs new isochrones at %.1f Gyr  (%d MC sets)'%(AGE_MYR/1000,n),
             fontsize=14,color=INK,x=0.005,ha='left',y=0.999)
fig.tight_layout(rect=[0.012,0.012,1,0.975])
out=f'{A}/fig6_og_vs_new_47tuc_{int(AGE_MYR/1000)}gyr.png'
fig.savefig(out,bbox_inches='tight'); print('wrote',os.path.basename(out))
