import glob, numpy as np, pandas as pd
d = pd.concat([pd.read_csv(f) for f in glob.glob('/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/verify-ast-correction-1/res_s*.csv')])
d['ks2'] = d.ks**2
order = ['shared','mirror','split','tq','hq','qsplit']
print('n per design', d.groupby('design').size().to_dict())
t = d.groupby(['split','design']).ks.mean().unstack()[order]; print((t*1e3).round(2).to_string())
m2 = d.groupby(['split','design']).ks2.mean().unstack()[order]
ex = m2.sub(m2['shared'], axis=0)   # AST-term excess in E[KS^2]
print('\nexcess E[KS^2] over shared (1e-5):'); print((ex*1e5).round(2).to_string())
U = {'mirror':1,'split':4,'tq':16/3,'hq':6,'qsplit':8}
pooled = ex.mean(); se = ex.std()/np.sqrt(len(ex))
print('\npooled excess (1e-5) and per-unit c:')
for k,u in U.items(): print(f'  {k:7s} units {u:5.2f} excess {pooled[k]*1e5:6.2f} +- {se[k]*1e5:5.2f}  per-unit {pooled[k]/u*1e5:6.2f}')
# fit excess = c*units across disjoint designs (through origin)
dj = ['split','tq','hq','qsplit']; uu = np.array([U[k] for k in dj]); yy = pooled[dj].values
c = (uu@yy)/(uu@uu); print(f'disjoint fit c = {c*1e5:.2f}e-5 per unit; predicted 1-unit excess {c*1e5:.2f}e-5; mirror measured {pooled["mirror"]*1e5:.2f}e-5')
# ratio of qsplit to split excess (predicted 2), per split
r = ex['qsplit']/ex['split']; print('qsplit/split excess ratio per split:', r.round(2).to_dict(), ' pooled', round(pooled['qsplit']/pooled['split'],2))
# predicted mirrored null from disjoint fit
sh2 = m2['shared'].mean(); print(f'sqrt(E shared^2)={np.sqrt(sh2):.5f}; predicted 1-unit null rms = {np.sqrt(sh2+c):.5f}; split rms {np.sqrt(sh2+4*c):.5f}')
pm = d.groupby('design').ks.agg(['mean','std']); pm['sig3']=pm['mean']+3*pm['std']; print(pm.loc[order].round(5).to_string())
