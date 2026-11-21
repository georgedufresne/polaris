import os, sys, tarfile, numpy as np, pandas as pd
sys.path.insert(0,'/dartfs-hpc/scratch/f007058/audit_ngc6791')
from trackio import read_track

T='/dartfs-hpc/scratch/f007058/Tracks_test'
WORK='/dartfs-hpc/scratch/f007058/audit_ngc6791/extracted'
VARS=[int(x) for x in open('/scratch/f007058/claude-550484/-dartfs-hpc/84db2d29-6b4d-4b69-a72e-63390c149db0/scratchpad/good_vars.txt').read().split()]
os.makedirs(WORK, exist_ok=True)

rows=[]
for v in VARS:
    d=os.path.join(WORK,str(v)); os.makedirs(d,exist_ok=True)
    for pref,tar in (('full',f'{T}/mcfeh031.{v}.tar'),('low',f'{T}/mcfehLOW031.{v}.tar')):
        if not os.path.exists(tar): continue
        with tarfile.open(tar) as tf: tf.extractall(d)
    for f in sorted(os.listdir(d)):
        if not f.endswith('.gz') or f.startswith('iso_'): continue
        pref = 'low' if f.startswith('ml') else 'full'
        try: mass=int(f.split('.')[0].lstrip('ml'))/1000.0
        except ValueError: continue
        t=read_track(os.path.join(d,f))
        if t is None:
            rows.append(dict(var=v,pass_=pref,mass=mass,n=0,status='UNREADABLE')); continue
        age=t['star_age']; lT=t['log_Teff']; lL=t['log_L']
        bad=(~np.isfinite(lT)).sum()+(~np.isfinite(lL)).sum()+(~np.isfinite(age)).sum()
        mono=bool(np.all(np.diff(age)>=0))
        st='OK'
        if bad: st='NONFINITE'
        elif not mono: st='AGE_NONMONOTONIC'
        elif len(age)<50: st='TOO_SHORT'
        elif lT.max()>5.5 or lT.min()<3.0: st='LOGT_OUT_OF_RANGE'
        elif lL.max()>6.0 or lL.min()<-4.0: st='LOGL_OUT_OF_RANGE'
        rows.append(dict(var=v,pass_=pref,mass=mass,n=len(age),
                         age_end=age[-1],logT_min=lT.min(),logT_max=lT.max(),
                         logL_min=lL.min(),logL_max=lL.max(),
                         nonfinite=bad,monotonic=mono,status=st))
df=pd.DataFrame(rows)
df.to_csv('/dartfs-hpc/scratch/f007058/audit_ngc6791/track_audit.csv',index=False)
print('tracks audited:',len(df))
print(df.status.value_counts().to_string())
print('\nper-pass counts:'); print(df.groupby('pass_').size().to_string())
