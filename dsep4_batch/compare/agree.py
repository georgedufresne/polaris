#!/usr/bin/env python3
"""Is dsepX_fixed vs dsep4 agreement within 1%?

Both grids use identical MC parameters and endage=1.8D10, so the only
difference is the code. Tracks are compared AT THE SAME AGE (the natural
common coordinate), over the age range both cover, and the fractional
differences in Teff and L are reported.

1% corresponds to 0.00432 dex in log Teff or log L.
"""
import gzip, os, sys, tarfile
import numpy as np, pandas as pd

XCOLS = ['star_age','star_mass','star_mdot','he_core_mass','c_core_mass',
         'log_L','log_LH','log_LHe','log_Teff','log_R','log_g',
         'surface_h1','surface_he3','surface_he4','surface_c12','surface_o16',
         'log_center_T','log_center_Rho','center_gamma',
         'center_h1','center_he4','center_c12']
PCT1 = np.log10(1.01)          # 0.004321 dex == 1%

def read_dsepX(buf):
    txt = buf.decode('utf-8', 'replace')
    lines = txt.split('\n')
    i = next((k+1 for k,l in enumerate(lines) if 'star_age' in l and 'log_Teff' in l), None)
    if i is None: return None
    v = []
    for l in lines[i:]:
        s = l.strip()
        if not s: continue
        try: v.append(float(s))
        except ValueError: continue
    n = len(v)//len(XCOLS)
    if n < 10: return None
    a = np.array(v[:n*len(XCOLS)]).reshape(n, len(XCOLS))
    return {c: a[:,j] for j,c in enumerate(XCOLS)}

def read_dsep4(p):
    with gzip.open(p,'rt',errors='replace') as fh: lines = fh.read().split('\n')
    i = next((k for k,l in enumerate(lines) if 'star_age' in l and 'log_Teff' in l), None)
    if i is None: return None
    h = lines[i].split(); rows=[]
    for l in lines[i+1:]:
        p2 = l.split()
        if len(p2)==len(h):
            try: rows.append([float(x) for x in p2])
            except ValueError: pass
    if len(rows) < 10: return None
    a = np.array(rows)
    return {c: a[:,j] for j,c in enumerate(h)}

rows=[]
for run in range(20000, 20011):
    tar = f'/dartfs-hpc/scratch/f007058/Tracks_noco/mcfeh031.{run}.tar'
    d4d = f'/dartfs-hpc/scratch/f007058/Tracks_dsep4_parts_18Gyr/{run}'
    if not (os.path.exists(tar) and os.path.isdir(d4d)): continue
    xm={}
    with tarfile.open(tar) as tf:
        for m in tf.getmembers():
            b=os.path.basename(m.name)
            if b.startswith('m') and b.endswith('.gz'):
                xm[b]=gzip.decompress(tf.extractfile(m).read())
    for key, buf in sorted(xm.items()):
        p4 = f'{d4d}/{key}'
        if not os.path.exists(p4): continue
        tx, t4 = read_dsepX(buf), read_dsep4(p4)
        if tx is None or t4 is None: continue
        mass = int(key[1:5])/1000.0
        # common age range, log-spaced so the long MS doesn't swamp the giant branch
        lo = max(tx['star_age'][1], t4['star_age'][1])
        hi = min(tx['star_age'][-1], t4['star_age'][-1])
        if not (hi > lo*1.01): continue
        ages = np.logspace(np.log10(lo), np.log10(hi), 400)
        f = lambda t,k: np.interp(ages, t['star_age'], t[k])
        dT = f(t4,'log_Teff') - f(tx,'log_Teff')
        dL = f(t4,'log_L')    - f(tx,'log_L')
        rows.append({'run':run,'mass':mass,
                     'medT':np.median(np.abs(dT)), 'medL':np.median(np.abs(dL)),
                     'fracT_1pct':np.mean(np.abs(dT)<PCT1),
                     'fracL_1pct':np.mean(np.abs(dL)<PCT1),
                     'biasT':np.median(dT), 'biasL':np.median(dL)})
df=pd.DataFrame(rows)
df.to_csv('/dartfs-hpc/scratch/f007058/dsep4_batch/compare/agreement.csv', index=False)
print(f"compared {len(df)} (run, mass) track pairs across "
      f"{df.run.nunique()} runs, {df.mass.nunique()} masses\n")
print(f"1% == {PCT1:.5f} dex\n")
print("ALL MASSES")
print(f"  median |dlogTeff| = {df.medT.median():.5f} dex  ({(10**df.medT.median()-1)*100:.2f}%)")
print(f"  median |dlogL|    = {df.medL.median():.5f} dex  ({(10**df.medL.median()-1)*100:.2f}%)")
print(f"  fraction of sampled points within 1% in Teff: {df.fracT_1pct.mean()*100:.1f}%")
print(f"  fraction within 1% in L                     : {df.fracL_1pct.mean()*100:.1f}%")
print(f"  pairs whose MEDIAN agrees within 1% in Teff : {(df.medT<PCT1).mean()*100:.1f}%")
print(f"  pairs whose MEDIAN agrees within 1% in L    : {(df.medL<PCT1).mean()*100:.1f}%")
print("\nBY MASS")
print(f"{'mass':>6}{'n':>4}{'med|dTeff|%':>13}{'med|dL|%':>11}{'%pts<1% T':>11}{'%pts<1% L':>11}")
for m,g in df.groupby('mass'):
    print(f"{m:>6.2f}{len(g):>4}{(10**g.medT.median()-1)*100:>13.2f}"
          f"{(10**g.medL.median()-1)*100:>11.2f}{g.fracT_1pct.mean()*100:>11.1f}{g.fracL_1pct.mean()*100:>11.1f}")
