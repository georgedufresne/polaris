#!/usr/bin/env python3
"""Three-way sCMD comparison: OG vs New vs Emulated (NGC 6791)."""
import numpy as np, pandas as pd
from scipy.stats import gaussian_kde

LEGS = {
 'OG':       '/dartfs-hpc/scratch/f007058/analysis_OG',
 'New':      '/dartfs-hpc/scratch/f007058/analysis_noco',
 'Emulated': '/dartfs-hpc/scratch/f007058/analysis_emul',
}
FILT = {'OG':'/dartfs-hpc/rc/home/8/f007058/filtered_ks_scmd.csv',
        'New':'/dartfs-hpc/scratch/f007058/analysis_noco/filtered_ks_scmd.csv',
        'Emulated':'/dartfs-hpc/scratch/f007058/analysis_emul/filtered_ks_scmd.csv'}
OUT='/dartfs-hpc/scratch/f007058/scmd_3way_comparison.txt'
N_BOOT=2000; RNG=np.random.default_rng(20261016)

boot=pd.read_csv('/dartfs-hpc/rc/home/8/f007058/Bootstrap/bootstrap_2d_ks_results.csv')
kde=gaussian_kde(boot.ks_statistic.values,bw_method='scott')
age_bins=np.arange(3.95,12.15,0.1); age_ctr=(age_bins[:-1]+age_bins[1:])/2

def load(p):
    d=pd.read_csv(p); d['w']=kde(d.ks_statistic.values)
    return d.loc[d.groupby(['simulation_id','age_gyr'])['ks_statistic'].idxmin()].reset_index(drop=True)

def wmean(v,w,bins,ctr):
    h,_=np.histogram(v,bins=bins,weights=w)
    if h.sum()<=0: return np.nan
    p=h/h.sum(); return float(np.sum(ctr*p))

S={k:pd.read_csv(f'{v}/scmd_only_parameter_summary.csv').set_index('Parameter') for k,v in LEGS.items()}
D={k:load(v) for k,v in FILT.items()}

L=[]; A=L.append
A('='*100)
A('NGC 6791 sCMD-ONLY THREE-WAY:  OG  vs  New (fixed physics)  vs  Emulated (CONF1DENCE, central z=0)')
A('='*100); A('')
A('Identical pipeline, identical KS cut (0.0468, 3 sigma of the original bootstrap), identical')
A('DM/E(B-V) grid. Only the isochrones differ.'); A('')
A(f'{"leg":<11}{"rows kept":>11}{"MC sets":>10}{"best KS":>10}')
for k in LEGS:
    d=D[k]; A(f'{k:<11}{len(d):>11,}{d.simulation_id.nunique():>10,}{d.ks_statistic.min():>10.4f}')
A('')
A('-'*100); A('PARAMETER POSTERIORS'); A('-'*100)
A(f'{"Parameter":<16}{"OG":>22}{"New":>22}{"Emulated":>22}')
for p in S['OG'].index:
    f=4 if p=='He (Y)' else 3
    row=f'{p:<16}'
    for k in LEGS:
        if p in S[k].index:
            r=S[k].loc[p]; row+=f'{f"{r.Mean:.{f}f} +/- {r.Std:.{f}f}":>22}'
        else: row+=f'{"-":>22}'
    A(row)
A('')
A('-'*100); A('AGE: paired difference on MC sets common to each pair'); A('-'*100)
pairs=[('New','OG'),('Emulated','OG'),('Emulated','New')]
A(f'{"comparison":<20}{"delta age":>12}{"SE":>9}{"z":>8}{"95% CI":>24}')
verdict=[]
for a,b in pairs:
    da,db=D[a],D[b]
    common=np.intersect1d(da.simulation_id.unique(),db.simulation_id.unique())
    ca=da[da.simulation_id.isin(common)]; cb=db[db.simulation_id.isin(common)]
    ga=ca.groupby('simulation_id').indices; gb=cb.groupby('simulation_id').indices
    va,wa=ca.age_gyr.values,ca.w.values; vb,wb=cb.age_gyr.values,cb.w.values
    obs=wmean(va,wa,age_bins,age_ctr)-wmean(vb,wb,age_bins,age_ctr)
    sets=np.asarray(sorted(common)); diffs=np.empty(N_BOOT)
    for i in range(N_BOOT):
        pick=RNG.integers(0,len(sets),len(sets))
        ia=np.concatenate([ga[sets[j]] for j in pick]); ib=np.concatenate([gb[sets[j]] for j in pick])
        diffs[i]=wmean(va[ia],wa[ia],age_bins,age_ctr)-wmean(vb[ib],wb[ib],age_bins,age_ctr)
    se=diffs.std(ddof=1); lo,hi=np.percentile(diffs,[2.5,97.5]); z=obs/se if se>0 else np.nan
    A(f'{f"{a} - {b}":<20}{obs:>+12.4f}{se:>9.4f}{z:>8.2f}{f"[{lo:+.4f}, {hi:+.4f}]":>24}   n={len(common)}')
    verdict.append((f'{a} - {b}',z))
A('')
A('-'*100); A('VERDICT (age)'); A('-'*100)
for n,z in verdict:
    A(f'  {n:<20} z = {z:+5.2f}   {"SIGNIFICANT" if abs(z)>=2 else "consistent with zero"}')
A(''); A('='*100)
t='\n'.join(L); open(OUT,'w').write(t+'\n'); print(t); print(f'\nSaved -> {OUT}')
