"""Bootstrap rows with one code tree: Gaia fixture (default flags) and 47 Tuc with shared ASTs.

usage: python boot_cmp.py CODE_ROOT OUT_JSON
"""
import argparse
import dataclasses
import inspect
import json
import os
import sys
import time

root, out = sys.argv[1], sys.argv[2]
sys.path.insert(0, root)

import numpy as np  # noqa: E402

from isochronetoolbox.scmd import Context, Setup  # noqa: E402
from isochronetoolbox.scmd import __main__ as M  # noqa: E402
from isochronetoolbox.scmd import bootstrap as B  # noqa: E402
from isochronetoolbox.scmd.pipeline import context  # noqa: E402

R = os.path.dirname(os.path.abspath(__file__))
params = inspect.signature(B.realisations).parameters
res = {"params": list(params)}

# Gaia (NGC 6791 fixture), default flags
fixture = os.path.join(root, "tests", "scmd", "data", "iso_20001_8.3Gyr.txt")
g = Context(Setup(iso=fixture, dm_grid=(13.27, 13.315, 0.01), ebv_grid=(0.12, 0.165, 0.01)))
t = time.time()
rows = B.realisations(g, 20001, 8.3, 13.28, 0.13, range(4))
res["gaia_default"] = [{k: (float(v) if isinstance(v, (float, np.floating)) else int(v)) for k, v in r.items()} for r in rows]
print(f"gaia default {time.time() - t:.1f} s", [r["ks_statistic"] for r in rows], flush=True)

# 47 Tuc from the old grid's setup.json, narrowed to the one cell
ns = argparse.Namespace(out=os.path.join(R, "olddir_nopreset"), preset=None, iso=None, obs=None, bc=None,
                        exact_age=False, chunk=None, ridge_width=None, ks_percentile=None)
s = dataclasses.replace(M._setup(ns), dm_grid=(13.38, 13.385, 0.01), ebv_grid=(0.015, 0.0175, 0.005))
ctx = context(s)
shared = {"split_ast": False} if "split_ast" in params else {"resample_ast": False}
t = time.time()
rows = B.realisations(ctx, 10099, 11.4, 13.38, 0.015, [0, 1], **shared)
res["tuc_shared"] = [{k: (float(v) if isinstance(v, (float, np.floating)) else int(v)) for k, v in r.items()} for r in rows]
print(f"47tuc shared {shared} {time.time() - t:.1f} s", [r["ks_statistic"] for r in rows], flush=True)
json.dump(res, open(out, "w"), indent=1)
