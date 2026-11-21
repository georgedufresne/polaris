import os,sys,numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
sys.path.insert(0,'/dartfs-hpc/scratch/f007058/audit_ngc6791')
from trackio import read_mist_iso

A='/dartfs-hpc/scratch/f007058/audit_ngc6791'
OG='/dartfs-hpc/scratch/f007058/subset_isochrones_ngc6791'
NEW=f'{A}/pyeep_work/isos'
VARS=[int(x) for x in open('/scratch/f007058/claude-550484/-dartfs-hpc/84db2d29-6b4d-4b69-a72e-63390c149db0/scratchpad/good_vars.txt').read().split()]
VARS=[v for v in VARS if os.path.exists(f'{NEW}/{v}/isochrones_{v}.txt')]
C_OG,C_NEW='#2a78d6','#eb6834'
plt.rcParams.update({'font.size':11,'axes.grid':True,'grid.alpha':0.25,'grid.linewidth':0.6,
 'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'#fcfcfb',
 'axes.facecolor':'#fcfcfb','savefig.dpi':150})
INK='#0b0b0b'; INK2='#52514e'

def grab(path, age):
    cols,a,_=read_mist_iso(path)
    if a is None: return None
    ag=a[:,cols.index('isochrone_age_yr')]
    sel=np.isclose(ag,age,rtol=1e-6)
    if sel.sum()==0: return None
    s=a[sel]
    return dict(eep=s[:,cols.index('EEP')], lT=s[:,cols.index('log_Teff')],
                lL=s[:,cols.index('log_L')], m=s[:,cols.index('initial_mass')])

res={}
fig,ax=plt.subplots(2,2,figsize=(13,10))
for j,age in enumerate([8.0e9,12.0e9]):
    dT=[];dL=[]
    for i,v in enumerate(VARS):
        o=grab(f'{OG}/isochrones_{v}.txt',age); n=grab(f'{NEW}/{v}/isochrones_{v}.txt',age)
        if o is None or n is None: continue
        lab_o='OG isochrones' if i==0 else None
        lab_n='New isochrones' if i==0 else None
        ax[0,j].plot(o['lT'],o['lL'],lw=1.1,color=C_OG,alpha=0.55,label=lab_o)
        ax[0,j].plot(n['lT'],n['lL'],lw=1.1,color=C_NEW,alpha=0.55,ls='--',label=lab_n)
        # residuals on the shared EEP grid
        common=np.intersect1d(o['eep'],n['eep'])
        if len(common)>5:
            oi={int(e):k for k,e in enumerate(o['eep'])}; ni={int(e):k for k,e in enumerate(n['eep'])}
            dT.append(np.array([n['lT'][ni[int(e)]]-o['lT'][oi[int(e)]] for e in common]))
            dL.append(np.array([n['lL'][ni[int(e)]]-o['lL'][oi[int(e)]] for e in common]))
            ax[1,j].plot(common,dT[-1],lw=0.9,color=C_NEW,alpha=0.45)
    ax[0,j].invert_xaxis()
    ax[0,j].set_xlabel('log $T_{\\rm eff}$'); ax[0,j].set_ylabel('log $L/L_\\odot$')
    ax[0,j].set_title(f'{"AB"[j]}. {age/1e9:.0f} Gyr isochrones — {len(VARS)} MC sets',loc='left',color=INK)
    ax[0,j].legend(frameon=False,loc='upper left')
    ax[1,j].axhline(0,color=INK2,lw=0.8)
    ax[1,j].set_xlabel('EEP'); ax[1,j].set_ylabel('$\\Delta$ log $T_{\\rm eff}$ (new $-$ OG)')
    ax[1,j].set_title(f'{"CD"[j]}. {age/1e9:.0f} Gyr residual vs OG',loc='left',color=INK)
    if dT:
        allT=np.concatenate(dT); allL=np.concatenate(dL)
        res[age]=(allT,allL)
        ax[1,j].text(0.02,0.04,'median $\\Delta$logT = %+.4f dex\nRMS = %.4f dex\n(%.0f K at 5500 K)'%(
            np.median(allT),np.sqrt(np.mean(allT**2)),np.sqrt(np.mean(allT**2))*np.log(10)*5500),
            transform=ax[1,j].transAxes,fontsize=10,color=INK2,va='bottom')
for a_ in ax.ravel(): a_.tick_params(colors=INK2)
fig.suptitle('New DSEP isochrones vs original — NGC 6791',fontsize=14,color=INK,x=0.01,ha='left')
fig.tight_layout(rect=[0,0,1,0.97])
fig.savefig(f'{A}/fig2_iso_new_vs_og.png',bbox_inches='tight')
print('wrote fig2_iso_new_vs_og.png  (vars used: %d)'%len(VARS))
for age,(t,l) in res.items():
    print('  %4.0f Gyr: dlogT med %+.4f rms %.4f | dlogL med %+.4f rms %.4f'%(
        age/1e9,np.median(t),np.sqrt(np.mean(t**2)),np.median(l),np.sqrt(np.mean(l**2))))
