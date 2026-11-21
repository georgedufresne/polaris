import sys,glob,os,numpy as np,pandas as pd; sys.path.insert(0,'.')
from trackio import read_mist_iso
rows=[]
for p in sorted(glob.glob('pyeep_work/isos/noco*/isochrones_noco*.txt')):
    r=int(p.split('noco')[-1].replace('.txt',''))
    sz=os.path.getsize(p)
    cols,a,_=read_mist_iso(p)
    if a is None: rows.append(dict(run=r,size=sz,rows=0,eep_max=0,logL_max=np.nan,ok=False)); continue
    e=a[:,cols.index('EEP')]; L=a[:,cols.index('log_L')]
    rows.append(dict(run=r,size=sz,rows=len(a),eep_max=int(e.max()),
                     logL_max=float(L.max()),ok=bool(e.max()>=600 and L.max()>2.0)))
d=pd.DataFrame(rows); d.to_csv('iso_audit_997.csv',index=False)
print('isochrones audited:',len(d))
print('  FULL  (EEP>=600 and logL>2):',d.ok.sum())
print('  BAD   (truncated)          :',(~d.ok).sum())
if (~d.ok).any():
    b=d[~d.ok].sort_values('run')
    print('\nbad runs:',' '.join(map(str,b.run.tolist())))
    print('\n  %-8s %10s %7s %8s %9s'%('run','size','rows','eep_max','logL_max'))
    for _,x in b.iterrows():
        print('  %-8d %10d %7d %8d %9.2f'%(x['run'],x['size'],x['rows'],x['eep_max'],x['logL_max']))
