"""Generate 47 Tuc sCMDs of one cell with one code tree (repo or dev).

usage: python gen.py CODE_ROOT OUT_NPZ VARIANT [--grid]

The Setup is built exactly as the CLI builds it (``__main__._setup``) from a
copy of the old grid's setup.json, then the variant's fields are replaced.
"""
import argparse
import dataclasses
import os
import sys
import time

root, out, variant = sys.argv[1], sys.argv[2], sys.argv[3]
do_grid = "--grid" in sys.argv
sys.path.insert(0, root)

import numpy as np  # noqa: E402

from isochronetoolbox.scmd import __main__ as M  # noqa: E402
from isochronetoolbox.scmd.pipeline import context, regenerate  # noqa: E402

R = os.path.dirname(os.path.abspath(__file__))
ns = argparse.Namespace(out=os.path.join(R, "olddir_nopreset"), preset=None, iso=None, obs=None, bc=None,
                        exact_age=False, chunk=None, ridge_width=None, ks_percentile=None)
s = M._setup(ns)
over = {"default": {}, "independent": {"ast_errors": "independent"}, "isochrone": {"cut": "isochrone"},
        "indep_iso": {"ast_errors": "independent", "cut": "isochrone"}}[variant]
s = dataclasses.replace(s, **over)
print(root, variant, "ast_errors", s.ast_errors, "cut", s.cut, {k: getattr(s, k, None) for k in
      ("ast_k", "ast_metric", "radial", "crn")}, flush=True)
t = time.time()
ctx = context(s)
print(f"context {time.time() - t:.1f} s", flush=True)
RUN, AGE = 10099, 11.4
idm = int(np.argmin(np.abs(s.dm - 13.38)))
ieb = int(np.argmin(np.abs(s.ebv - 0.015)))
gdm, geb = float(s.dm[idm]), float(s.ebv[ieb])
print("grid dm", repr(gdm), "grid ebv", repr(geb), flush=True)
_, make = ctx.model(RUN, AGE)
res = {}
for seed in (0, 1, 4242):
    t = time.time()
    cat = make(np.random.RandomState(seed), gdm, geb)
    res[f"seed{seed}"] = cat
    ks = ctx.stage3(cat)["ks_statistic"]
    res[f"ks_seed{seed}"] = np.array(ks)
    print(f"seed {seed}: {cat.shape} KS {ks!r}  {time.time() - t:.1f} s", flush=True)
# internals
res["dist"] = ctx.ast.dist
res["cdf"] = ctx.ast.cdf
res["nbr"] = ctx.ast.nearest.nbr
res["n_target"] = np.array(ctx.n_target)
if "--stream" in sys.argv:
    import itertools
    t = time.time()
    for j, (d, e, cat) in enumerate(itertools.islice(ctx.grid(RUN, AGE)[1], 3)):
        res[f"stream{j}"] = cat
        res[f"stream{j}_de"] = np.array([d, e])
    print(f"stream: 3 cells {time.time() - t:.1f} s", flush=True)
if do_grid:
    t = time.time()
    cat = regenerate(ctx, RUN, AGE, 13.38, 0.015)
    res["regen"] = cat
    ks = ctx.stage3(cat)["ks_statistic"]
    res["ks_regen"] = np.array(ks)
    print(f"regenerate: {cat.shape} KS {ks!r}  {time.time() - t:.1f} s", flush=True)
np.savez(out, **res)
print("saved", out, flush=True)
