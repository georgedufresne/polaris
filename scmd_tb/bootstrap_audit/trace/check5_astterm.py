"""How large is the AST-list term?  KS between two fresh sCMDs of the best cell whose errors come from
different AST lists.  argv[1] = config: shared (full/full), split (A/B), half_full (A/full), full_half (full/B);
argv[2] = number of pairs.  Linear model: KS^2 ~ s^2 + u * (list-difference variance in units of the full list):
shared 0, split 4, half_full 1, full_half 1; data vs model 1."""
import sys, json, time, numpy as np
sys.path.insert(0, '/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox')
from isochronetoolbox.scmd import Setup, context
from isochronetoolbox.scmd.pipeline import _window
from isochronetoolbox.scmd.stats import ks_2d, rank2d
cfg, npairs = sys.argv[1], int(sys.argv[2])
G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
S = Setup(**{k: v for k, v in json.load(open(G + '/setup.json')).items() if k in Setup.__dataclass_fields__})
ctx = context(S)
full = ctx.ast
A, B = full.halves()
lists = {"shared": (full, full), "split": (A, B), "half_full": (A, full), "full_half": (full, B)}[cfg]
_, make = ctx.model(10099, 11.4)
def scmd(ast, sd):
    ctx.ast = ast
    try:
        p = ctx.points(make(np.random.RandomState(sd), 13.38, 0.015))
    finally:
        ctx.ast = full
    return p[_window(p, S.ks_window)]
ks = []
for k in range(npairs):
    ref = scmd(lists[0], 5000 + k)
    sim = scmd(lists[1], 7000 + k)
    v = ks_2d(ref, sim, rank2d(ref) / len(ref))[0]
    ks.append(v)
    print(f"{cfg} pair {k}: KS {v:.5f}", flush=True)
ks = np.array(ks)
print(f"RESULT {cfg}: mean {ks.mean():.5f} sd {ks.std(ddof=1):.5f} sem {ks.std(ddof=1)/np.sqrt(len(ks)):.5f} n {len(ks)}")
