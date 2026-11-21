import os,io,gzip,tarfile,sys,numpy as np,pandas as pd
sys.path.insert(0,'/dartfs-hpc/scratch/f007058/audit_ngc6791')
from trackio import TRACK_COLS
T='/dartfs-hpc/scratch/f007058/Tracks_tuc_noco'
V='/dartfs-hpc/scratch/f007058/47Tuc_var'
def rd(fh):
    L=fh.read().decode('utf-8','replace').split('\n'); st=None
    for i,l in enumerate(L):
        if 'star_age' in l and 'log_Teff' in l: st=i+1; break
    if st is None: return None
    v=[]
    for l in L[st:]:
        s=l.strip()
        if s:
            try: v.append(float(s))
            except ValueError: pass
    n=len(v)//22
    return None if n==0 else {c:np.array(v[:n*22]).reshape(n,22)[:,j] for j,c in enumerate(TRACK_COLS)}
rows=[]
for p in sorted(os.listdir(T)):
    if not p.startswith('mcfeh077.'): continue
    run=int(p.split('.')[1])
    d={}
    vf=f'{V}/varfeh-077.{run}'
    if os.path.exists(vf):
        for line in open(vf):
            q=line.split()
            if len(q)>=2: d[q[1]]=float(q[0])
    with tarfile.open(f'{T}/{p}') as tf:
        for mem in tf.getmembers():
            nm=os.path.basename(mem.name)
            if not nm.startswith('m') or nm.startswith('iso_') or not nm.endswith('.gz'): continue
            try: mass=int(nm.split('.')[0].lstrip('m'))/1000.
            except ValueError: continue
            t=rd(gzip.GzipFile(fileobj=io.BytesIO(tf.extractfile(mem).read())))
            if t is None: rows.append(dict(run=run,mass=mass,status='UNREADABLE')); continue
            lL=t['log_L']; a=t['star_age']
            rows.append(dict(run=run,mass=mass,n=len(a),end_gyr=a[-1]/1e9,
                logL_max=lL.max(),capped=len(a)>=59990,reach_tip=bool(lL.max()>2.3),
                aMLT=d.get('CMIXLA',np.nan),FeH=d.get('FeH',np.nan),
                status='OK' if np.isfinite(lL).all() else 'NONFINITE'))
df=pd.DataFrame(rows); df.to_csv('/dartfs-hpc/scratch/f007058/audit_ngc6791/tuc_track_audit.csv',index=False)
print('tracks:',len(df)); print(df.status.value_counts().to_string())
ok=df[df.status.eq('OK')]
print('capped (>=60k): %d'%ok.capped.sum())
print('reach RGB tip : %d (%.1f%%)'%(ok.reach_tip.sum(),100*ok.reach_tip.mean()))
print('max models used: %d'%ok.n.max())
