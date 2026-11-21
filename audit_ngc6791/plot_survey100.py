import os,pickle,numpy as np,pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

A='/dartfs-hpc/scratch/f007058/audit_ngc6791'
df=pd.read_csv(f'{A}/survey100_tracks.csv')
n_bad=df.status.ne('OK').sum()
df=df[df.status.eq('OK')].copy()   # drop UNREADABLE (NaN params break corr)
print('dropped %d unreadable tracks'%n_bad)
curves=pickle.load(open(f'{A}/survey100_curves.pkl','rb'))
VARS=sorted(curves.keys())
plt.rcParams.update({'font.size':8,'axes.grid':True,'grid.alpha':0.2,'grid.linewidth':0.4,
 'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'#fcfcfb',
 'axes.facecolor':'#fcfcfb','savefig.dpi':110})
INK='#0b0b0b'; INK2='#52514e'
OKC='#1baf7a'; BADC='#eb6834'

# ---------- FIG 3: grid of HR diagrams, one panel per MC set ----------
n=len(VARS); ncol=10; nrow=int(np.ceil(n/ncol))
fig,axes=plt.subplots(nrow,ncol,figsize=(ncol*2.1,nrow*2.0),sharex=True,sharey=True)
axes=np.atleast_2d(axes)
cmap=plt.get_cmap('viridis')
for i,v in enumerate(VARS):
    ax=axes[i//ncol,i%ncol]
    s=df[df['var']==v]
    cap=int(np.sum(s.capped.values.astype(bool))); aC=s.aC.iloc[0]
    for mass,lT,lL in curves[v]:
        ax.plot(lT,lL,lw=0.4,color=cmap(min(1,(mass-0.2)/2.8)),alpha=0.85)
    edge = BADC if cap>=5 else (OKC if cap==0 else INK2)
    for sp in ax.spines.values(): sp.set_color(edge); sp.set_linewidth(1.6)
    ax.spines['top'].set_visible(True); ax.spines['right'].set_visible(True)
    ax.set_title('%d  $\\alpha_C$=%.3f\ncapped %d/39'%(v,aC,cap),fontsize=6.5,
                 color=(BADC if cap>=5 else INK),pad=2)
    ax.tick_params(labelsize=5,colors=INK2)
for j in range(n,nrow*ncol): axes[j//ncol,j%ncol].axis('off')
axes[0,0].invert_xaxis()
fig.supxlabel('log $T_{\\rm eff}$',color=INK); fig.supylabel('log $L/L_\\odot$',color=INK)
leg=[Line2D([],[],color=OKC,lw=2,label='0 capped tracks'),
     Line2D([],[],color=INK2,lw=2,label='1-4 capped'),
     Line2D([],[],color=BADC,lw=2,label='>=5 capped (RGB at risk)')]
fig.legend(handles=leg,loc='upper right',frameon=False,fontsize=9,ncol=3,
           bbox_to_anchor=(0.995,0.998))
fig.suptitle('NGC 6791 — evolutionary tracks, first %d MC sets (39 masses each; panel border = health)'%n,
             fontsize=13,color=INK,x=0.005,ha='left',y=0.999)
fig.tight_layout(rect=[0.012,0.012,1,0.982])
fig.savefig(f'{A}/fig3_tracks_100sets.png',bbox_inches='tight'); plt.close(fig)
print('wrote fig3_tracks_100sets.png')

# ---------- FIG 4: the diagnosis ----------
g=df.groupby('var').agg(aC=('aC','first'),aMLT=('aMLT','first'),FeH=('FeH','first'),
                        capped=('capped','sum'),tip=('reach_tip','sum')).reset_index()
# groupby on bool-from-CSV yields object dtype, which breaks .corr()
for c in ('aC','aMLT','FeH','capped','tip'): g[c]=g[c].astype(float)
fig,ax=plt.subplots(1,3,figsize=(15,4.4))
plt.rcParams.update({'font.size':11})
ax[0].scatter(g.aC,g.capped,s=26,color=BADC,alpha=0.8,edgecolor='none')
ax[0].set_xlabel('core overshoot $\\alpha_C$'); ax[0].set_ylabel('capped tracks per set (of 39)')
ax[0].set_title('A. nmod cap is driven by core overshoot\n r = %+.3f'%g.aC.corr(g.capped),
                loc='left',color=INK)
ax[1].scatter(g.aC,g.tip,s=26,color=OKC,alpha=0.8,edgecolor='none')
ax[1].set_xlabel('core overshoot $\\alpha_C$'); ax[1].set_ylabel('tracks reaching RGB tip (of 39)')
ax[1].set_title('B. and it removes the RGB\n r = %+.3f'%g.aC.corr(g.tip),loc='left',color=INK)
ax[2].scatter(g.aMLT,g.capped,s=26,color=INK2,alpha=0.7,edgecolor='none')
ax[2].set_xlabel('mixing length $\\alpha_{MLT}$'); ax[2].set_ylabel('capped tracks per set')
ax[2].set_title('C. control: no $\\alpha_{MLT}$ dependence\n r = %+.3f'%g.aMLT.corr(g.capped),
                loc='left',color=INK)
for a in ax: a.tick_params(colors=INK2); a.grid(alpha=0.25,lw=0.6)
fig.suptitle('Why some isochrones lose their RGB — %d MC sets'%len(g),fontsize=13,color=INK,x=0.005,ha='left')
fig.tight_layout(rect=[0,0,1,0.93])
fig.savefig(f'{A}/fig4_cap_diagnosis.png',bbox_inches='tight')
print('wrote fig4_cap_diagnosis.png')
print('\nsets with >=5 capped: %d/%d (%.0f%%)'%((g.capped>=5).sum(),len(g),100*(g.capped>=5).mean()))
print('sets with 0 capped : %d/%d'%((g.capped==0).sum(),len(g)))
