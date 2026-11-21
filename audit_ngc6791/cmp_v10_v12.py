#!/usr/bin/env python3
"""Compare pyeep v1.0 (stock input.nml) against v1.2 (input_dsee.nml) isochrones
built from the same NGC 6791 tracks."""
import numpy as np, pandas as pd, os, warnings
warnings.filterwarnings('ignore')
A='/dartfs-hpc/scratch/f007058/audit_ngc6791'
W=f'{A}/pyeep_work'
ld=lambda p: pd.read_csv(p,sep=r'\s+',comment='#',header=None,usecols=[0,1,2,7,10],
                         names=['EEP','age','m0','logL','logTeff'],engine='c').dropna()
AGE=8.5e9
# the sCMD fit window, in log L (M_G 3.00-4.72 at DM 13.28, BC_G=-0.15)
WLO,WHI=0.073,0.761
runs=[l.strip() for l in open(f'{A}/v12_runs.txt') if l.strip()]
rows=[]
for R in runs:
    p10=f'{W}/isos/noco{R}/isochrones_noco{R}.txt'
    p12=f'{W}/isos_v12/v12{R}/isochrones_v12{R}.txt'
    if not (os.path.exists(p10) and os.path.getsize(p12 if os.path.exists(p12) else p10)>0): continue
    if not os.path.exists(p12) or os.path.getsize(p12)==0: continue
    r={'run':R}
    for tag,p in (('v10',p10),('v12',p12)):
        d=ld(p); ag=np.unique(d['age']); a=ag[np.argmin(abs(ag-AGE))]
        s=d[d['age']==a]
        r[f'n_{tag}']=len(s); r[f'eeplo_{tag}']=int(s.EEP.min()); r[f'eephi_{tag}']=int(s.EEP.max())
        r[f'Llo_{tag}']=s.logL.min(); r[f'Lhi_{tag}']=s.logL.max()
        r[f'nages_{tag}']=len(ag)
        r[f'cov_{tag}']=(s.logL.max()>=WHI) and (s.logL.min()<=WLO)
        r[f'size_{tag}']=os.path.getsize(p)/1e6
    rows.append(r)
df=pd.DataFrame(rows)
print(f"N runs compared: {len(df)}   (8.5 Gyr isochrone from each)\n")
print(f"{'':<22}{'v1.0':>22}{'v1.2':>22}")
def line(lab,a,b,f='{:.2f}'):
    print(f"{lab:<22}{f.format(a):>22}{f.format(b):>22}")
line('ages per file',df.nages_v10.mean(),df.nages_v12.mean(),'{:.0f}')
line('points at 8.5 Gyr',df.n_v10.mean(),df.n_v12.mean(),'{:.0f}')
line('first EEP',df.eeplo_v10.mean(),df.eeplo_v12.mean(),'{:.0f}')
line('last EEP',df.eephi_v10.mean(),df.eephi_v12.mean(),'{:.0f}')
line('min log L',df.Llo_v10.mean(),df.Llo_v12.mean(),'{:+.3f}')
line('max log L',df.Lhi_v10.mean(),df.Lhi_v12.mean(),'{:+.3f}')
line('file size (MB)',df.size_v10.mean(),df.size_v12.mean(),'{:.1f}')
print()
print(f"last-EEP spread : v1.0 {df.eephi_v10.min()}-{df.eephi_v10.max()}   "
      f"v1.2 {df.eephi_v12.min()}-{df.eephi_v12.max()}")
print(f"max log L spread: v1.0 {df.Lhi_v10.std():.4f}            v1.2 {df.Lhi_v12.std():.4f}")
print()
print(f"sCMD fit window is log L {WLO:.3f} .. {WHI:.3f}")
print(f"  isochrones spanning the FULL window:  v1.0 {df.cov_v10.sum()}/{len(df)}   "
      f"v1.2 {df.cov_v12.sum()}/{len(df)}")
short=df[~df.cov_v12]
if len(short):
    print(f"  v1.2 falls short of log L {WHI:.3f} by a median of "
          f"{(WHI-short.Lhi_v12).median():.3f} dex")
print()
print("Per-run detail (first 12):")
print(f"{'run':>7}{'v1.0 EEP':>12}{'v1.2 EEP':>12}{'v1.0 maxL':>11}{'v1.2 maxL':>11}")
for _,r in df.head(12).iterrows():
    print(f"{r.run:>7}{str(r.eeplo_v10)+'-'+str(r.eephi_v10):>12}"
          f"{str(r.eeplo_v12)+'-'+str(r.eephi_v12):>12}{r.Lhi_v10:>+11.3f}{r.Lhi_v12:>+11.3f}")
df.to_csv(f'{A}/v10_vs_v12_summary.csv',index=False)
print(f"\nSaved -> {A}/v10_vs_v12_summary.csv")
