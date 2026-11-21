import os,sys,numpy as np,pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
sys.path.insert(0,'/dartfs-hpc/scratch/f007058/audit_ngc6791')
from trackio import read_track

A='/dartfs-hpc/scratch/f007058/audit_ngc6791'
VARS=[int(x) for x in open('/scratch/f007058/claude-550484/-dartfs-hpc/84db2d29-6b4d-4b69-a72e-63390c149db0/scratchpad/good_vars.txt').read().split()]
df=pd.read_csv(f'{A}/track_audit.csv')

plt.rcParams.update({'font.size':11,'axes.grid':True,'grid.alpha':0.25,
    'grid.linewidth':0.6,'axes.spines.top':False,'axes.spines.right':False,
    'figure.facecolor':'#fcfcfb','axes.facecolor':'#fcfcfb','savefig.dpi':150})
INK='#0b0b0b'; INK2='#52514e'

fig,ax=plt.subplots(2,2,figsize=(13,10))

# --- A: HR diagram, all masses, one var set (sequential ramp by mass) ---
v0=VARS[0]; d=f'{A}/extracted/{v0}'
files=sorted([f for f in os.listdir(d) if f.endswith('.gz') and not f.startswith('iso_') and f.startswith('m') and not f.startswith('ml')])
cmap=plt.get_cmap('viridis')
masses=[]
for f in files:
    m=int(f.split('.')[0].lstrip('m'))/1000.0; masses.append(m)
mn,mx=min(masses),max(masses)
for f,m in zip(files,masses):
    t=read_track(os.path.join(d,f))
    if t is None: continue
    ax[0,0].plot(t['log_Teff'],t['log_L'],lw=0.9,color=cmap((m-mn)/(mx-mn)),alpha=0.9)
ax[0,0].invert_xaxis()
ax[0,0].set_xlabel('log $T_{\\rm eff}$'); ax[0,0].set_ylabel('log $L/L_\\odot$')
ax[0,0].set_title(f'A. Evolutionary tracks, var {v0} (39 masses)',loc='left',color=INK)
sm=plt.cm.ScalarMappable(cmap=cmap,norm=plt.Normalize(mn,mx))
cb=fig.colorbar(sm,ax=ax[0,0]); cb.set_label('initial mass ($M_\\odot$)',color=INK2)

# --- B: models per track vs mass ---
full=df[df.pass_=='full']
for v in VARS:
    s=full[full['var']==v].sort_values('mass')
    ax[0,1].plot(s.mass,s.n,lw=1.0,color='#2a78d6',alpha=0.35)
ax[0,1].set_xscale('log'); ax[0,1].set_xlabel('initial mass ($M_\\odot$)')
ax[0,1].set_ylabel('models per track')
ax[0,1].set_title('B. Track length vs mass (11 sets, full pass)',loc='left',color=INK)

# --- C: final age vs mass ---
for v in VARS:
    s=full[full['var']==v].sort_values('mass')
    ax[1,0].plot(s.mass,s.age_end/1e9,lw=1.0,color='#eb6834',alpha=0.35)
ax[1,0].set_xscale('log'); ax[1,0].set_yscale('log')
ax[1,0].set_xlabel('initial mass ($M_\\odot$)'); ax[1,0].set_ylabel('final age (Gyr)')
ax[1,0].set_title('C. Track end age vs mass',loc='left',color=INK)

# --- D: RGB tip luminosity vs mass ---
for v in VARS:
    s=full[full['var']==v].sort_values('mass')
    ax[1,1].plot(s.mass,s.logL_max,lw=1.0,color='#1baf7a',alpha=0.35)
ax[1,1].set_xscale('log'); ax[1,1].set_xlabel('initial mass ($M_\\odot$)')
ax[1,1].set_ylabel('max log $L/L_\\odot$ (tip)')
ax[1,1].set_title('D. Peak luminosity vs mass',loc='left',color=INK)

for a in ax.ravel(): a.tick_params(colors=INK2)
fig.suptitle(f'New DSEP tracks — audit overview  ({len(df)} tracks, {df.status.eq("OK").sum()} OK)',
             fontsize=14,color=INK,x=0.01,ha='left')
fig.tight_layout(rect=[0,0,1,0.97])
fig.savefig(f'{A}/fig1_track_audit.png',bbox_inches='tight')
print('wrote fig1_track_audit.png')
