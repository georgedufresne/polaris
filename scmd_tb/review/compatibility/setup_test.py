"""setup.json compatibility through the CLI's _setup.

usage: python setup_test.py CODE_ROOT
"""
import argparse
import dataclasses
import json
import os
import shutil
import sys

root = sys.argv[1]
sys.path.insert(0, root)
from isochronetoolbox.scmd import __main__ as M  # noqa: E402

R = os.path.dirname(os.path.abspath(__file__))
tag = os.path.basename(root)


def ns(out, preset=None):
    return argparse.Namespace(out=out, preset=preset, iso=None, obs=None, bc=None, exact_age=False, chunk=None,
                              ridge_width=None, ks_percentile=None)


# 1. old grid directory (setup.json without the new fields), with and without --preset 47tuc
for preset in (None, "47tuc"):
    d = os.path.join(R, f"olddir_{'preset' if preset else 'nopreset'}")
    try:
        s = M._setup(ns(d, preset))
        print(f"[{tag}] old dir, preset={preset}: OK", {k: getattr(s, k, "<no field>") for k in
              ("ast_errors", "ast_k", "ast_metric", "radial", "crn", "crn_seed", "ks_percentile")})
    except SystemExit as e:
        print(f"[{tag}] old dir, preset={preset}: EXIT {e}")

# 2. a new directory from --preset 47tuc: which keys get written
d = os.path.join(R, f"newdir_{tag}")
shutil.rmtree(d, ignore_errors=True)
s = M._setup(ns(d, "47tuc"))
saved = json.load(open(os.path.join(d, "setup.json")))
print(f"[{tag}] new dir keys added vs old grid's setup.json:",
      sorted(set(saved) - set(json.load(open(os.path.join(R, "olddir_nopreset", "setup.json"))))))

# 3. a setup.json written for the new model (smooth / crowding / corrected / crn), read by this tree
d = os.path.join(R, f"newmodel_setup_{tag}")
os.makedirs(d, exist_ok=True)
shutil.copy("/dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/newmodel/setup.json", d)
s = M._setup(ns(d))
print(f"[{tag}] new-model setup.json read as:", {k: getattr(s, k, "<dropped>") for k in
      ("ast_errors", "ast_k", "ast_metric", "radial", "crn", "crn_seed")})
try:
    M._setup(ns(d, "47tuc"))
    print(f"[{tag}] new-model setup.json with --preset 47tuc: OK")
except SystemExit as e:
    print(f"[{tag}] new-model setup.json with --preset 47tuc: EXIT {str(e)[-160:]}")
