"""The grid's own generator (code snapshot tuc_grid_nearest_bin02/code, copied to ./gridcode):
1. regenerate the best cell exactly as the grid drew it (grid seed, cell 257 of the
   sequence) and score it against the data with ks_2d -> should give the grid's 0.016828;
2. save its KS-window points (reference G0);
3. data vs N fresh realisations of the same cell with the grid generator
   (np.random.RandomState(j), j = 0..N-1, full AST list).

usage: gridcode_check.py N
"""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/reference-selection")
from common import AGE, EVI, MU, OUT, RUN, Log, load_ctx, prep_target, setup_path, ks_many, window_points  # noqa: E402

setup_path("gridcode")
import isochronetoolbox  # noqa: E402
from isochronetoolbox.scmd.ast import make_scmd_ast  # noqa: E402
from isochronetoolbox.scmd.model import apparent, select_age  # noqa: E402
from isochronetoolbox.scmd.pipeline import regenerate  # noqa: E402
from isochronetoolbox.scmd.stats import ks_2d  # noqa: E402

N = int(sys.argv[1])
log = Log("gridcode_check")
log("toolbox imported from", isochronetoolbox.__file__)
S, ctx = load_ctx()
log("context ready")

cat = regenerate(ctx, RUN, AGE, MU, EVI)
g0 = window_points(ctx, S, cat)
np.save(f"{OUT}/ref_G0.npy", g0)
ks_g0 = ks_2d(ctx.obs3, g0, ctx.obs3_cdf)[0]
grid = pd.read_csv("/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02/combined_2d_ks_results.csv")
row = grid[(grid.simulation_id == RUN) & np.isclose(grid.age_gyr, AGE) & np.isclose(grid.distance_modulus, MU)
           & np.isclose(grid.extinction_ebv, EVI)]
log(f"G0 = grid's sCMD regenerated with the grid code: n {len(g0)}, KS vs data {ks_g0:.7f}; "
    f"grid table {row.ks_statistic.values}")

# the grid's generator for this cell (ASTContext.grid of the snapshot, without the shared seed)
age_yr, rows = select_age(ctx.isochrones(RUN), AGE, S.exact_age)
absmag = ctx.absolute(rows)
mass = rows["initial_mass"].values
coef = {"F606W": 1.0, "F814W": 0.0}


def make(rs, dm, ebv):
    return make_scmd_ast(rs, mass, apparent(absmag, dm, ebv, coef, "F606W"), ctx.ast, ctx.fid,
                         ctx.n_target, ctx.v_window, S.pad, S.imf_slope, S.binary_fraction,
                         S.martin_sign, S.cut, S.fid_width, rcut=ctx.rcut, errors=S.ast_errors)


data = prep_target(ctx.obs3, ctx.obs3_cdf)
out = []
for j in range(N):
    p = window_points(ctx, S, make(np.random.RandomState(j), MU, EVI))
    r = ks_many(p, [data])[0]
    out.append({"seed": j, "n": len(p), **{f"data_{k}": v for k, v in r.items()}})
    log(f"grid-code realisation {j}: KS vs data {r['ks']:.6f}")
df = pd.DataFrame(out)
df.to_csv(f"{OUT}/gridcode_data_vs_fresh.csv", index=False)
log(f"data vs {N} grid-code realisations: mean {df.data_ks.mean():.6f} std {df.data_ks.std(ddof=1):.6f}")
log("done")
