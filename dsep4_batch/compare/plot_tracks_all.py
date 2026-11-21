#!/usr/bin/env python3
"""dsepX_fixed vs dsep4 tracks, runs 20000-20010. Same MC params, same 18 Gyr endage."""
import gzip, os, tarfile
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

XC = ['star_age','star_mass','star_mdot','he_core_mass','c_core_mass','log_L','log_LH',
      'log_LHe','log_Teff','log_R','log_g','surface_h1','surface_he3','surface_he4',
      'surface_c12','surface_o16','log_center_T','log_center_Rho','center_gamma',
      'center_h1','center_he4','center_c12']
OUT='/dartfs-hpc/scratch/f007058/dsep4_batch/compare'

def rx(buf):
    lines=buf.decode('utf-8','replace').split('\n')
    i=next((k+1 for k,l in enumerate(lines) if 'star_age' in l and 'log_Teff' in l),None)
    if i is None: return None
    v=[]
    for l in lines[i:]:
        s=l.strip()
        if not s: continue
        try: v.append(float(s))
        except ValueError: continue
    n=len(v)//len(XC)
    if n<10: return None
    a=np.array(v[:n*len(XC)]).reshape(n,len(XC))
    return {c:a[:,j] for j,c in enumerate(XC)}

def r4(p):
    with gzip.open(p,'rt',errors='replace') as fh: lines=fh.read().split('\n')
    i=next((k for k,l in enumerate(lines) if 'star_age' in l and 'log_Teff' in l),None)
    if i is None: return None
    h=lines[i].split(); rows=[]
    for l in lines[i+1:]:
        p2=l.split()
        if len(p2)==len(h):
            try: rows.append([float(x) for x in p2])
            except ValueError: pass
    if len(rows)<10: return None
    a=np.array(rows); return {c:a[:,j] for j,c in enumerate(h)}

RUNS=[str(r) for r in range(20000,20011)]
# ---- Fig 1: every mass, one panel per run ----
fig,axes=plt.subplots(3,4,figsize=(22,15))
for ax,run in zip(axes.ravel(),RUNS):
    tar=f'/dartfs-hpc/scratch/f007058/Tracks_noco/mcfeh031.{run}.tar'
    d4=f'/dartfs-hpc/scratch/f007058/Tracks_dsep4_parts_18Gyr/{run}'
    xm={}
    if os.path.exists(tar):
        with tarfile.open(tar) as tf:
            for m in tf.getmembers():
                b=os.path.basename(m.name)
                if b.startswith('m') and b.endswith('.gz'):
                    xm[b]=gzip.decompress(tf.extractfile(m).read())
    nx=n4=0
    for key in sorted(xm):
        t=rx(xm[key])
        if t is not None:
            ax.plot(t['log_Teff'],t['log_L'],'-',color='#1f77b4',lw=.7,alpha=.85); nx+=1
    if os.path.isdir(d4):
        for f in sorted(os.listdir(d4)):
            if not (f.startswith('m') and f.endswith('.gz')): continue
            t=r4(f'{d4}/{f}')
            if t is not None:
                ax.plot(t['log_Teff'],t['log_L'],'-',color='#d62728',lw=.7,alpha=.85); n4+=1
    ax.set_title(f'run {run}   dsepX {nx} / dsep4 {n4}',fontsize=10)
    ax.invert_xaxis(); ax.grid(alpha=.25)
    ax.set_xlabel('log $T_{eff}$',fontsize=8); ax.set_ylabel('log $L/L_\\odot$',fontsize=8)
    ax.tick_params(labelsize=8)
for ax in axes.ravel()[len(RUNS):]: ax.axis('off')
from matplotlib.lines import Line2D
axes.ravel()[0].legend(handles=[Line2D([],[],color='#1f77b4',label='dsepX_fixed'),
                                Line2D([],[],color='#d62728',label='dsep4')],fontsize=9)
fig.suptitle('NGC 6791 tracks, runs 20000-20010 — dsepX_fixed vs dsep4 '
             '(identical MC parameters, endage 18 Gyr)',fontsize=16)
fig.tight_layout(rect=[0,0,1,0.965]); fig.savefig(f'{OUT}/fig_tracks_all_runs.png',dpi=100); plt.close(fig)
print('saved fig_tracks_all_runs.png')
