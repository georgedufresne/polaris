import sys, json, time, numpy as np, pandas as pd
sys.path.insert(0, '/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox')
from isochronetoolbox.scmd import Setup, context
from isochronetoolbox.scmd.pipeline import _window
from isochronetoolbox.scmd.stats import ks_2d, rank2d
G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
split_seed = int(sys.argv[1]); refseeds = [int(x) for x in sys.argv[2].split(',')]; K = int(sys.argv[3]); tag = sys.argv[4]
S = Setup(**{k: v for k, v in json.load(open(G + '/setup.json')).items() if k in Setup.__dataclass_fields__})
ctx = context(S); full = ctx.ast
A, B = full.halves(seed=split_seed)
_, make = ctx.model(10099, 11.4)
def cat(lst, seed):
    ctx.ast = lst
    try:
        p = ctx.points(make(np.random.RandomState(seed), 13.38, 0.015))
    finally:
        ctx.ast = full
    return p[_window(p, S.ks_window)]
t0 = time.time()
reals = {'full': [cat(full, i) for i in range(K)], 'B': [cat(B, i) for i in range(K)]}
print('reals', time.time() - t0, [len(r) for r in reals['full'][:3]], flush=True)
rows = []
designs = {'shared': ('full', 'full'), 'split': ('A', 'B'), 'mirror': ('A', 'full'), 'rmirror': ('full', 'B')}
for rs in refseeds:
    refs = {'full': cat(full, rs), 'A': cat(A, rs)}
    cdfs = {k: rank2d(v) / len(v) for k, v in refs.items()}
    for d, (rl, ml) in designs.items():
        for i in range(K):
            ks = ks_2d(refs[rl], reals[ml][i], cdfs[rl])[0]
            rows.append(dict(split_seed=split_seed, ref_seed=rs, design=d, real_seed=i, ks=ks, nref=len(refs[rl]), nsim=len(reals[ml][i])))
        print(rs, d, np.mean([r['ks'] for r in rows[-K:]]), time.time() - t0, flush=True)
    pd.DataFrame(rows).to_csv(f'/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/verify-ast-correction-0/res_{tag}.csv', index=False)
