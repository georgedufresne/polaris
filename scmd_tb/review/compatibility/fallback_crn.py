"""(a) repo code + new-model setup.json: is ast_errors="smooth" silently the independent model?
(b) dev code + crn: does regenerate() return the grid cell's catalogue bit for bit?

usage: python fallback_crn.py CODE_ROOT {fallback|crn}
"""
import argparse
import dataclasses
import os
import sys
import time

root, what = sys.argv[1], sys.argv[2]
sys.path.insert(0, root)
import numpy as np  # noqa: E402

from isochronetoolbox.scmd import __main__ as M  # noqa: E402
from isochronetoolbox.scmd.pipeline import context, regenerate  # noqa: E402

R = os.path.dirname(os.path.abspath(__file__))
tag = os.path.basename(root)
ns = argparse.Namespace(out=os.path.join(R, f"newmodel_setup_{tag}"), preset=None, iso=None, obs=None, bc=None,
                        exact_age=False, chunk=None, ridge_width=None, ks_percentile=None)
s = M._setup(ns)
RUN, AGE = 10099, 11.4
gdm = float(s.dm[int(np.argmin(np.abs(s.dm - 13.38)))])
geb = float(s.ebv[int(np.argmin(np.abs(s.ebv - 0.015)))])

if what == "fallback":
    out = {}
    for errs in ("smooth", "independent"):
        ctx = context(dataclasses.replace(s, ast_errors=errs))
        _, make = ctx.model(RUN, AGE)
        out[errs] = make(np.random.RandomState(0), gdm, geb)
    print(f"[{tag}] ast_errors='smooth' catalogue == ast_errors='independent' catalogue:",
          np.array_equal(out["smooth"], out["independent"]))
else:
    print(f"[{tag}] crn={s.crn} crn_seed={s.crn_seed} ast_errors={s.ast_errors} metric={s.ast_metric} radial={s.radial}")
    ctx = context(s)
    _, make = ctx.model(RUN, AGE)
    t = time.time()
    grid_cell = make(np.random.RandomState(s.crn_seed), gdm, geb)      # what ctx.grid yields for the cell
    # confirm against ctx.grid itself (crn: each cell is independent, so islice to it)
    import itertools
    k = int(np.argmin(np.abs(s.dm - 13.38))) * len(s.ebv) + int(np.argmin(np.abs(s.ebv - 0.015)))
    d, e, from_grid = next(itertools.islice(ctx.grid(RUN, AGE)[1], k, k + 1))
    print("grid() cell dm, ebv:", repr(d), repr(e), " equals make(gdm, geb):", np.array_equal(from_grid, grid_cell))
    rg = regenerate(ctx, RUN, AGE, 13.38, 0.015)                       # the CSV / cell.json values
    print(f"regenerate(13.38, 0.015) bitwise equal to the grid cell: {np.array_equal(rg, grid_cell)}; "
          f"max |diff| {np.max(np.abs(rg - grid_cell)) if rg.shape == grid_cell.shape else 'shape differs'}; "
          f"shapes {rg.shape} {grid_cell.shape}")
    a, b = ctx.stage3(grid_cell)["ks_statistic"], ctx.stage3(rg)["ks_statistic"]
    print(f"KS grid cell {a!r}  KS regenerate {b!r}  equal: {a == b}   ({time.time() - t:.0f} s)")
