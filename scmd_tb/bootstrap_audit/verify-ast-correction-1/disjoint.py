"""Disjoint-list scaling test: does the AST-list term in E[KS^2] scale as 1/n_ref + 1/n_model
(the premise behind 'split = 4 units, data-vs-model = 1 unit')?  Uses only disjoint lists, plus mirror for comparison."""
import sys, json, time, numpy as np, pandas as pd
sys.path.insert(0, '/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox')
from isochronetoolbox.scmd import Setup, context
from isochronetoolbox.scmd.pipeline import _window
from isochronetoolbox.scmd.stats import ks_2d, rank2d
from isochronetoolbox.scmd.ast import ArtificialStars
from scipy.interpolate import LinearNDInterpolator
G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
OUT = '/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/verify-ast-correction-1'
s, rseed, K = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
T0 = time.time()
def log(m): print(f'[s{s} {time.time()-T0:6.1f}s] {m}', flush=True)
S = Setup(**{k: v for k, v in json.load(open(G + '/setup.json')).items() if k in Setup.__dataclass_fields__})
ctx = context(S); full = ctx.ast
_, make = ctx.model(10099, 11.4)
n = len(pd.read_csv(full._args[0], usecols=['InputX']))
u = np.random.RandomState(s).random_sample(n)
L = {'F': full}
masks = {'HA': u < 0.5, 'HB': u >= 0.5, 'Q1': u < 0.25, 'Q3': (u >= 0.5) & (u < 0.75), 'T3': ~((u >= 0.5) & (u < 0.75))}
for k, m in masks.items():
    L[k] = ArtificialStars(*full._args, rows=m)
hA, hB = full.halves(s)
assert hA.n_ast == L['HA'].n_ast and np.array_equal(hA.nearest.nbr, L['HA'].nearest.nbr)
log('lists ' + ' '.join(f'{k}:{v.n_ast}' for k, v in L.items()))
def cat(lst, seed):
    ctx.ast = lst
    try:
        p = ctx.points(make(np.random.RandomState(seed), 13.38, 0.015))
    finally:
        ctx.ast = full
    return p[_window(p, S.ks_window)]
refs = {k: cat(L[k], rseed) for k in ('F', 'HA', 'Q1', 'T3')}
cdfs = {k: rank2d(v) / len(v) for k, v in refs.items()}
assert ctx.ast is full
log('refs ' + str({k: len(v) for k, v in refs.items()}))
# design: name -> (ref list, realisation list, units of 1/N_full in Var(AST term) under 1/N scaling)
D = {'shared': ('F', 'F', 0.0), 'mirror': ('HA', 'F', 1.0), 'split': ('HA', 'HB', 4.0),
     'tq': ('T3', 'Q3', 4/3 + 4), 'hq': ('HA', 'Q3', 6.0), 'qsplit': ('Q1', 'Q3', 8.0)}
rows = []
for i in range(K):
    for rl in ('F', 'HB', 'Q3'):
        sim = cat(L[rl], i)
        f = LinearNDInterpolator(sim, rank2d(sim) / len(sim), fill_value=0)
        for d, (a, b, un) in D.items():
            if b != rl: continue
            ks = float(np.percentile(np.abs(cdfs[a] - f(refs[a])), 95))
            if i == 0 and d in ('split', 'qsplit'):
                chk = ks_2d(refs[a], sim, cdfs[a])[0]
                log(f'check vs ks_2d {d}: {ks:.10f} {chk:.10f}')
            rows.append(dict(split=s, ref_seed=rseed, design=d, units=un, seed=i, ks=ks, n_ref=len(refs[a]), n_sim=len(sim)))
    pd.DataFrame(rows).to_csv(f'{OUT}/res_s{s}_r{rseed}.csv', index=False)
    log(f'i={i} ' + ' '.join(f"{r['design']} {r['ks']:.5f}" for r in rows[-6:]))
log('done')
