"""Data-vs-model KS of fresh realisations of the best cell under one code version.
argv[1] = code root (latest repo or the grid's frozen copy), argv[2] = label."""
import sys, json, time, numpy as np
root, label = sys.argv[1], sys.argv[2]
sys.path.insert(0, root)
import isochronetoolbox.scmd as M
from isochronetoolbox.scmd import Setup, context
from isochronetoolbox.scmd.pipeline import _window, select_age
from isochronetoolbox.scmd.model import apparent
from isochronetoolbox.scmd.ast import make_scmd_ast
from isochronetoolbox.scmd.stats import ks_2d
print(label, "package at", M.__file__)
G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
S = Setup(**{k: v for k, v in json.load(open(G + '/setup.json')).items() if k in Setup.__dataclass_fields__})
ctx = context(S)
age_yr, rows = select_age(ctx.isochrones(10099), 11.4, S.exact_age)
absmag = ctx.absolute(rows); mass = rows["initial_mass"].values
coef = {"F606W": 1.0, "F814W": 0.0}
sizes = []
orig = ctx.rcut.trace
ctx.rcut.trace = lambda c, m: (sizes.append(len(c)), orig(c, m))[1]
out = []
for i in range(int(sys.argv[3]) if len(sys.argv) > 3 else 6):
    rs = np.random.RandomState(i)
    cat = make_scmd_ast(rs, mass, apparent(absmag, 13.38, 0.015, coef, "F606W"), ctx.ast, ctx.fid, ctx.n_target,
                        ctx.v_window, S.pad, S.imf_slope, S.binary_fraction, S.martin_sign, S.cut, S.fid_width,
                        rcut=ctx.rcut, errors=S.ast_errors)
    p = ctx.points(cat); p = p[_window(p, S.ks_window)]
    ks = ks_2d(ctx.obs3, p, ctx.obs3_cdf)[0]
    out.append(ks)
    print(f"{label} seed {i}: n {len(p)} trace-input {sizes[-1]} KS {ks:.5f}", flush=True)
out = np.array(out)
print(f"{label}: data-vs-model KS mean {out.mean():.5f} std {out.std(ddof=1):.5f} (n={len(out)})")
