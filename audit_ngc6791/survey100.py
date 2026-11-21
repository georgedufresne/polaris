"""Audit + plot evolutionary tracks for the first N complete MC sets."""
import os,sys,io,tarfile,gzip,numpy as np,pandas as pd
sys.path.insert(0,'/dartfs-hpc/scratch/f007058/audit_ngc6791')
from trackio import TRACK_COLS

T='/dartfs-hpc/scratch/f007058/Tracks_test'
A='/dartfs-hpc/scratch/f007058/audit_ngc6791'
VARD='/dartfs-hpc/rc/home/8/f007058/mcdsep/MC_var'
NWANT=int(sys.argv[1]) if len(sys.argv)>1 else 100

def read_stream(fh):
    lines=fh.read().decode('utf-8','replace').split('\n')
    st=None
    for i,l in enumerate(lines):
        if 'star_age' in l and 'log_Teff' in l: st=i+1; break
    if st is None: return None
    vals=[]
    for l in lines[st:]:
        s=l.strip()
        if not s: continue
        try: vals.append(float(s))
        except ValueError: pass
    n=len(vals)//22
    if n==0: return None
    a=np.array(vals[:n*22]).reshape(n,22)
    return {c:a[:,j] for j,c in enumerate(TRACK_COLS)}

# pick the first NWANT vars with a complete 39-member full-pass tar
vars_=[]
for v in range(20000,21200):
    p=f'{T}/mcfeh031.{v}.tar'
    if not os.path.exists(p): continue
    try:
        with tarfile.open(p) as tf:
            n=sum(1 for m in tf.getnames() if m.startswith('m') and not m.startswith('iso_'))
    except Exception: continue
    if n==39: vars_.append(v)
    if len(vars_)>=NWANT: break
print('sets selected:',len(vars_), vars_[0],'..',vars_[-1],flush=True)

rows=[]; curves={}
for k,v in enumerate(vars_):
    d={}
    for line in open(f'{VARD}/varfeh+031.{v}'):
        p=line.split()
        if len(p)>=2: d[p[1]]=float(p[0])
    cv=[]
    with tarfile.open(f'{T}/mcfeh031.{v}.tar') as tf:
        for mem in tf.getmembers():
            nm=os.path.basename(mem.name)
            if not nm.startswith('m') or nm.startswith('iso_') or not nm.endswith('.gz'): continue
            try: mass=int(nm.split('.')[0].lstrip('m'))/1000.
            except ValueError: continue
            t=read_stream(gzip.GzipFile(fileobj=io.BytesIO(tf.extractfile(mem).read())))
            if t is None:
                rows.append(dict(var=v,mass=mass,status='UNREADABLE')); continue
            age=t['star_age']; lL=t['log_L']; lT=t['log_Teff']; he=t['log_LHe']
            nmod=len(age)
            rows.append(dict(var=v,mass=mass,aC=d['ALPHAC'],aMLT=d['CMIXLA'],
                FeH=d['FeH'],KTTAU=int(d['KTTAU']),n=nmod,
                end_gyr=age[-1]/1e9,logL_max=lL.max(),logLHe_max=he.max(),
                capped=nmod>=9990, reach_tip=bool(lL.max()>2.3),
                reach13=bool(age[-1]>=13e9), he_ign=bool(he.max()>-3),
                status='OK' if np.isfinite(lL).all() and np.isfinite(lT).all() else 'NONFINITE'))
            cv.append((mass,lT[::6],lL[::6]))
    curves[v]=cv
    if (k+1)%10==0: print('  %d/%d sets'%(k+1,len(vars_)),flush=True)

df=pd.DataFrame(rows); df.to_csv(f'{A}/survey100_tracks.csv',index=False)
np.save(f'{A}/survey100_curves.npy',np.array([1]))  # marker
import pickle
pickle.dump(curves,open(f'{A}/survey100_curves.pkl','wb'))
print('\ntracks:',len(df))
print(df.status.value_counts().to_string())
print('capped: %d (%.1f%%)'%(df.capped.sum(),100*df.capped.mean()))
print('reach RGB tip: %d (%.1f%%)'%(df.reach_tip.sum(),100*df.reach_tip.mean()))
print('He ignited: %d'%df.he_ign.sum())
print('reach 13 Gyr: %d'%df.reach13.sum())
