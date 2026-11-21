"""Round 3: does the AST-centre re-estimate drive the split-half AST term?
Halves rebuilt with the FULL list's centre (NearestFixed), split seed 0 (the table split) and seed 2,
fresh refs (seeds 1000, 1005 for s=0; 1002 for s=2) -> compare with the same designs with per-half centres."""
import copy, os, sys
from common import *
from scipy.spatial import cKDTree
role, K = sys.argv[1], int(sys.argv[2])
out_csv = os.path.join(OUT, f'results3_{role}.csv')
T0 = time.time(); rows = []
def log(m): print(f'[{role} {time.time() - T0:7.1f}s] {m}', flush=True)
def save(design, tag, seed, ks, n_sim, n_ref):
    rows.append(dict(role=role, design=design, tag=str(tag), seed=int(seed), ks=ks, n_sim=n_sim, n_ref=n_ref))
    pd.DataFrame(rows).to_csv(out_csv, index=False)

class NearestFixed(NearestErrors):
    """NearestErrors with a given centre instead of one re-estimated from this list's ASTs."""
    def __init__(self, a, xc, yc, k=30, scale_mag=0.1, scale_r=200.0):
        x, y = a["InputX"].values, a["InputY"].values
        vin, iin = a["InputF606W"].values.astype(float), a["InputF814W"].values.astype(float)
        vout, iout = a["OutputF606W"].values.astype(float), a["OutputF814W"].values.astype(float)
        r = np.hypot(x - xc, y - yc)
        self.rec = ((np.hypot(x - a["OutputX"].values, y - a["OutputY"].values) < 0.5)
                    & (np.abs(vin - vout) < 0.75) & (np.abs(iin - iout) < 0.75))
        self.ev, self.ei = vin - vout, iin - iout
        self.k, self.scale = k, np.array([scale_mag, scale_mag, scale_r])
        pts = np.column_stack([vin, iin, r]) / self.scale
        self.tree = cKDTree(pts)
        self.nbr = self.tree.query(pts, k=k)[1].astype(np.int32)

S = setup(); ctx = context(S); full = ctx.ast; _, make = ctx.model(RUN, AGE)
a = pd.read_csv(S.ast, dtype=np.float32)
xc, yc = centre(a.InputX.values, a.InputY.values)
chk = NearestFixed(a, xc, yc)
log(f'full-list AST centre ({xc:.2f}, {yc:.2f}); NearestFixed(full) reproduces full.nearest: '
    f'{np.array_equal(chk.nbr, full.nearest.nbr) and np.array_equal(chk.rec, full.nearest.rec)}')
# observed fit stars' centre (used for the synthetic stars' radii) vs the AST centre
obs = ctx.obs
fx, fy = None, None
try:
    raw = pd.read_csv(S.obs_raw, sep=r'\s+')
    keep = ctx.rcut.select(raw['VIvega'].values, raw['Vvega'].values, ctx.rcut.ridge)
    fx, fy = centre(raw['x'].values[keep], raw['y'].values[keep])
    nfit = int(((raw['x'].values[keep] - 3000.0) ** 2 + (raw['y'].values[keep] - 3000.0) ** 2 <= 200.0 ** 2).sum())
    log(f'fit-star centre ({fx:.2f}, {fy:.2f}) from {nfit} fit stars within 200 px; offset from AST centre {np.hypot(fx - xc, fy - yc):.1f} px')
except Exception as e:
    log(f'fit-star centre not computed: {e}')

jobs = [(0, [1000, 1005]), (2, [1002])]
for s, seeds in jobs:
    h = np.random.RandomState(s).random_sample(len(a)) < 0.5
    hA, hB = copy.copy(full), copy.copy(full)
    hA.nearest = NearestFixed(a[h].reset_index(drop=True), xc, yc)
    hB.nearest = NearestFixed(a[~h].reset_index(drop=True), xc, yc)
    refs = [(f'{s}r{q}', Ref(sample(ctx, make, hA, q), f'A{s}fix_{q}')) for q in seeds]
    for i in range(K):
        sim = sample(ctx, make, hB, i)
        ks = ks_multi(sim, [r for _, r in refs])
        for (tag, r), k in zip(refs, ks):
            save('split_fixedcentre', tag, i, k, len(sim), len(r.pts))
        log(f's={s} i={i}: ' + ' '.join(f'{t} {k:.5f}' for (t, _), k in zip(refs, ks)))
log('done')
