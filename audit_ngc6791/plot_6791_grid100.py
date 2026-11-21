"""NGC 6791: OG vs new, BOTH ages in every panel. Colour = source, style = age."""
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
 'axes.facecolor':'#fcfcfb','savefig.dpi':110})
INK='#0b0b0b'; INK2='#52514e'; BAD='#b00020'

def at(path,age):
    cols,a,_=read_mist_iso(path)
    if a is None: return None
    ag=a[:,cols.index('isochrone_age_yr')]; s=a[np.isclose(ag,age,rtol=1e-6)]
    if len(s)==0: return None
    return s[:,cols.index('log_Teff')],s[:,cols.index('log_L')]

panels=[]
for r in runs:
    p=f'{A}/pyeep_work/isos/noco{r}/isochrones_noco{r}.txt'
    if not os.path.exists(p): continue
    d={}
    for A_,tag in [(8.0e9,'8'),(12.0e9,'12')]:
        d['og'+tag]=at(f'{OG}/isochrones_{r}.txt',A_)
        d['nw'+tag]=at(p,A_)
    v={}
    for line in open(f'/dartfs-hpc/rc/home/8/f007058/mcdsep/MC_var_noco/varfeh+031.{r}'):
        q=line.split()
        if len(q)>=2: v[q[1]]=float(q[0])
    # flag a gutted new isochrone (lost its giant branch)
    bad=all(d['nw'+t] is None or d['nw'+t][1].max()<1.0 for t in ('8','12'))
    panels.append((r,d,v,bad))

N=len(panels); ncol=10; nrow=int(np.ceil(N/ncol))
fig,axes=plt.subplots(nrow,ncol,figsize=(ncol*2.2,nrow*2.2),sharex=True,sharey=True)
axes=np.atleast_2d(axes)
nbad=0
for i,(r,d,v,bad) in enumerate(panels):
    ax=axes[i//ncol,i%ncol]
    for k,(c,ls) in {'og8':(C_OG,'-'),'og12':(C_OG,'--'),
                     'nw8':(C_NEW,'-'),'nw12':(C_NEW,'--')}.items():
        if d[k] is not None: ax.plot(d[k][0],d[k][1],lw=0.9,color=c,ls=ls,alpha=0.9)
    if bad:
        nbad+=1
        for sp in ax.spines.values(): sp.set_color(BAD); sp.set_linewidth(1.8)
        ax.spines['top'].set_visible(True); ax.spines['right'].set_visible(True)
    ax.set_title('%d  $\\alpha$=%.2f  %+.2f'%(r,v['CMIXLA'],v['FeH']),
                 fontsize=6.5,color=(BAD if bad else INK),pad=2)
    ax.tick_params(labelsize=5,colors=INK2)
for j in range(N,nrow*ncol): axes[j//ncol,j%ncol].axis('off')
axes[0,0].invert_xaxis()
fig.supxlabel('log $T_{\\rm eff}$',color=INK); fig.supylabel('log $L/L_\\odot$',color=INK)
fig.legend(handles=[Line2D([],[],color=C_OG,lw=2,label='OG'),
                    Line2D([],[],color=C_NEW,lw=2,label='New (overshoot off)'),
                    Line2D([],[],color=INK2,lw=2,ls='-',label='8 Gyr'),
                    Line2D([],[],color=INK2,lw=2,ls='--',label='12 Gyr'),
                    Line2D([],[],color=BAD,lw=2,label='new isochrone gutted')],
           loc='upper right',frameon=False,fontsize=10,ncol=5,bbox_to_anchor=(0.997,0.999))
fig.suptitle('NGC 6791 — OG vs new isochrones, 8 and 12 Gyr  (%d MC sets; %d gutted)'%(N,nbad),
             fontsize=14,color=INK,x=0.004,ha='left',y=0.999)
fig.tight_layout(rect=[0.010,0.010,1,0.978])
out=f'{A}/fig9_6791_grid100_both_ages.png'
fig.savefig(out,bbox_inches='tight')
print('wrote %s  (%d panels, %d gutted)'%(os.path.basename(out),N,nbad))
