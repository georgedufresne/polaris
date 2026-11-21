"""NGC 6791: CONF1DENCE emulator vs OG isochrones vs NEW (fixed-physics) isochrones.

Emulation, medoid sampling and physical closure follow
CONF1DENCE/scripts/emulator_compare/iso_compare_m15.py.
"""
import os, sys, json, argparse
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import torch, zuko

sys.path.insert(0, '/dartfs-hpc/rc/home/8/f007058/CONF1DENCE')
sys.path.insert(0, '/dartfs-hpc/scratch/f007058/audit_ngc6791')
import CONF1DENCE
from CONF1DENCE.utils.helper import data_transform as dt
from trackio import read_mist_iso, read_dsep_iso

A   = '/dartfs-hpc/scratch/f007058/audit_ngc6791'
OGD = '/dartfs-hpc/scratch/f007058/47Tuc_iso'
NEWD= f'{A}/pyeep_work/isos'
VARD= '/dartfs-hpc/scratch/f007058/47Tuc_var'
REPO= CONF1DENCE.__abspath__
CKPT= os.path.join(REPO,'data/model/DSEE_run2_lr1e4/best.pt')
FMT = os.path.join(REPO,'data/format/DSEE_iso_format.json')

LOG_G_SUN=4.4380; LOG_T_SUN=np.log10(5772.0)
_medoid_scales={'Log_G':0.06,'Log_T':0.02,'Log_L':0.06,'Log_R':0.06}
C_OG,C_NEW,C_EMU='#2a78d6','#eb6834','#1baf7a'

class Emu:
    def __init__(self, ckpt, fmt_path, device):
        self.device=device
        with open(fmt_path) as f: self.fmt=json.load(f)
        nf=sum(2 if s=='mantissa_exponent' else 1 for s in self.fmt['data_scaler'])
        nc=sum(2 if s=='mantissa_exponent' else 1 for s in self.fmt['label_scaler'])
        self.flow=zuko.flows.NSF(features=nf,context=nc,transforms=10,
                                 hidden_features=[256]*10).to(device)
        raw=torch.load(ckpt,map_location=device,weights_only=False)
        raw=raw['model'] if isinstance(raw,dict) and 'model' in raw else raw
        st={}
        for k,v in raw.items():
            if not torch.is_tensor(v): continue
            k=k[len('flow.'):] if k.startswith('flow.') else k
            if k.startswith(('transform.','base.')): st[k]=v
        miss,unexp=self.flow.load_state_dict(st,strict=False)
        bad=[k for k in miss if k.endswith(('weight','bias'))]
        if bad or unexp: raise RuntimeError(f'state_dict mismatch: {bad[:4]} {list(unexp)[:4]}')
        self.flow.eval()
        print(f'>> emulator loaded: features={nf} context={nc}')

def parse_varfeh(path):
    d={}
    for line in open(path):
        p=line.split()
        if len(p)>=2: d[p[1]]=float(p[0])
    return d

def physics_block(v):
    return np.array([v['FeH'], v.get('Yprim', v.get('He')), v['alphafe'],
        v['CMIXLA'], v['FGRY'], v['FGRZ'], v['KTTAU'], v['ALPHAE'], v['ALPHAC'],
        v['SSTANDARD(1)'],v['SSTANDARD(2)'],v['SSTANDARD(3)'],v['SSTANDARD(4)'],
        v['SSTANDARD(5)'],v['SSTANDARD(6)'],v['SSTANDARD(7)'], 1.0,
        v['alexcoef'],v['opalcoef2'],v['talphacoef'],v['plascoef'],v['cocoef']],
        dtype=np.float64)

def _sample_T(fd,n,temp):
    z=fd.base.rsample((n,))
    if temp!=1.0:
        mu=fd.base.base_dist.loc; z=mu+temp*(z-mu)
    return fd.transform.inv(z)

def representative_sample(s,names,scales=None):
    med=s.median(dim=0).values; diff=s-med.unsqueeze(0)
    if scales is not None:
        sv=torch.tensor([scales.get(n,1.0) for n in names],device=s.device,dtype=s.dtype)
        diff=diff/sv.view(1,1,-1)
    idx=(diff**2).sum(dim=2).argmin(dim=0)
    return s[idx,torch.arange(s.shape[1],device=s.device),:]

def emulate(emu,phys,age,args):
    fmt=emu.fmt
    m1=np.linspace(args.mass_low,args.mass_high,args.sample_size)
    labels=np.hstack((np.tile(phys,(args.sample_size,1)),
                      np.full((args.sample_size,1),age), m1.reshape(-1,1)))
    out=[]
    for i in range(0,len(labels),args.batch_labels):
        ch=labels[i:i+args.batch_labels].astype(np.float32)
        norm=dt.norm_data_torch(torch.tensor(ch),fmt['label_scaler'],fmt['label_scale'])
        norm=norm.to(emu.device)
        with torch.no_grad():
            fd=emu.flow(norm)
            sr=_sample_T(fd,args.n_samples,args.latent_temp)
            ns,nb,nfz=sr.shape
            ph=dt.inverse_norm_data_torch(sr.reshape(ns*nb,nfz),
                 fmt['data_scaler'],fmt['data_scale']).reshape(ns,nb,nfz)
            out.append(representative_sample(ph,fmt['data_header'],_medoid_scales).cpu().numpy())
    iso=pd.DataFrame(np.vstack(out),columns=fmt['data_header']); iso['mass']=m1
    logR=0.5*(iso['Log_L'].values-4.0*(iso['Log_T'].values-LOG_T_SUN))
    iso['Log_R']=logR; iso['Log_G']=LOG_G_SUN+np.log10(iso['mass'].values)-2.0*logR
    k=((iso.Log_T>3)&(iso.Log_T<4.5)&(iso.Log_L>-10)&(iso.Log_L<10))
    return iso[k].sort_values('mass').reset_index(drop=True)

def grab(path,age):
    """MIST-format (new, from pyeep) isochrone at a given age in yr."""
    cols,a,_=read_mist_iso(path)
    if a is None: return None
    ag=a[:,cols.index('isochrone_age_yr')]; s=a[np.isclose(ag,age,rtol=1e-6)]
    if len(s)==0: return None
    return pd.DataFrame({'Log_T':s[:,cols.index('log_Teff')],
                         'Log_L':s[:,cols.index('log_L')],
                         'mass':s[:,cols.index('initial_mass')]})

def grab_og(path,age):
    """DSEP-native OG isochrone; its header ages are in Myr."""
    d=read_dsep_iso(path)
    key=min(d.keys(),key=lambda a: abs(a-age/1e6)) if d else None
    if key is None or abs(key-age/1e6)>1.0: return None
    b=d[key]
    return pd.DataFrame({'Log_T':b['Log_T'],'Log_L':b['Log_L'],'mass':b['mass']})

def nn_rms(ref,test):
    """RMS nearest-neighbour distance in (logT,logL) from every REFERENCE point
    to the TEST curve.  Direction matters: measuring reference->test does not
    penalise `test` for extending beyond the reference (e.g. the emulator
    continues up the RGB past the real tip, which is expected behaviour)."""
    if ref is None or test is None or len(ref)==0 or len(test)==0: return np.nan
    pr=np.column_stack([ref.Log_T,ref.Log_L]); pt=np.column_stack([test.Log_T,test.Log_L])
    d=np.sqrt(((pr[:,None,:]-pt[None,:,:])**2).sum(-1)).min(axis=1)
    return float(np.sqrt(np.mean(d**2)))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--device',default=None)
    ap.add_argument('--sample-size',type=int,default=800)
    ap.add_argument('--n-samples',type=int,default=128)
    ap.add_argument('--mass-low',type=float,default=0.20)
    ap.add_argument('--mass-high',type=float,default=1.15)
    ap.add_argument('--latent-temp',type=float,default=1.0)
    ap.add_argument('--batch-labels',type=int,default=50_000)
    ap.add_argument('--vars',default='/dartfs-hpc/scratch/f007058/audit_ngc6791/tuc_good.txt')
    args=ap.parse_args()
    dev=args.device or ('cuda' if torch.cuda.is_available() else 'cpu')
    print('device:',dev)
    VARS=[int(x) for x in open(args.vars).read().split()]
    emu=Emu(CKPT,FMT,dev)
    pdir=f'{A}/plots_per_var_tuc'; os.makedirs(pdir,exist_ok=True)
    plt.rcParams.update({'font.size':11,'axes.grid':True,'grid.alpha':0.25,
      'grid.linewidth':0.6,'axes.spines.top':False,'axes.spines.right':False,
      'figure.facecolor':'#fcfcfb','axes.facecolor':'#fcfcfb','savefig.dpi':150})
    rows=[]
    for v in VARS:
        vf=parse_varfeh(f'{VARD}/varfeh-077.{v}')
        phys=physics_block(vf)
        fig,ax=plt.subplots(1,2,figsize=(13,5.4))
        for j,age in enumerate([6.0e9,8.6e9]):
            og=grab_og(f'{OGD}/feh077.{v}',age)
            nw=grab(f'{NEWD}/tuc{v}/isochrones_tuc{v}.txt',age)
            em=emulate(emu,phys,age,args)
            for d,c,ls,lab in ((og,C_OG,'-','OG isochrone'),
                               (nw,C_NEW,'--','New isochrone'),
                               (em,C_EMU,'-.','CONF1DENCE emulated')):
                if d is not None and len(d):
                    ax[j].plot(d.Log_T,d.Log_L,lw=1.6,color=c,ls=ls,label=lab,alpha=0.9)
            ax[j].invert_xaxis(); ax[j].set_xlabel('log $T_{\\rm eff}$')
            ax[j].set_ylabel('log $L/L_\\odot$')
            ax[j].set_title(f'{age/1e9:.0f} Gyr',loc='left',color='#0b0b0b')
            ax[j].legend(frameon=False,fontsize=9,loc='upper left')
            rows.append(dict(var=v,age_gyr=age/1e9,
                n_og=0 if og is None else len(og), n_new=0 if nw is None else len(nw),
                n_emu=len(em),
                rms_new_vs_og=nn_rms(og,nw), rms_emu_vs_og=nn_rms(og,em),
                rms_emu_vs_new=nn_rms(nw,em)))
        fig.suptitle(f'47 Tuc  MC var {v}   [Fe/H]={vf["FeH"]:+.3f}  '
                     f'Y={vf.get("Yprim",vf.get("He")):.3f}  '
                     f'$\\alpha_{{MLT}}$={vf["CMIXLA"]:.3f}  '
                     f'$\\alpha_{{C}}$={vf["ALPHAC"]:.3f}',
                     fontsize=13,color='#0b0b0b',x=0.01,ha='left')
        fig.tight_layout(rect=[0,0,1,0.93])
        fig.savefig(f'{pdir}/tuc_{v}_3way.png',bbox_inches='tight'); plt.close(fig)
        print(f'  var {v} -> tuc_{v}_3way.png')
    df=pd.DataFrame(rows); df.to_csv(f'{A}/compare_3way_metrics_tuc.csv',index=False)
    print('\n=== HRD nearest-neighbour RMS (dex) ===')
    print(df.to_string(index=False,float_format=lambda x:'%.4f'%x))

if __name__=='__main__': main()
