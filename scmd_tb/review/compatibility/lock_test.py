"""cell.json lock of the dev bootstrap against cell.json files written by older code.

usage: python lock_test.py CODE_ROOT
"""
import json
import os
import sys
import tempfile

root = sys.argv[1]
sys.path.insert(0, root)

from isochronetoolbox.scmd import Context, Setup  # noqa: E402
from isochronetoolbox.scmd import bootstrap as B  # noqa: E402

fixture = os.path.join(root, "tests", "scmd", "data", "iso_20001_8.3Gyr.txt")
ctx = Context(Setup(iso=fixture, dm_grid=(13.27, 13.315, 0.01), ebv_grid=(0.12, 0.165, 0.01)))
CELL = {"simulation_id": 20001, "age_gyr": 8.3, "distance_modulus": 13.28, "extinction_ebv": 0.13}
base = {**CELL, "n": 2}
cases = {
    "gaia, written by repo code (split_ast false)": {**base, "split_ast": False},
    "written by the oldest code (no flag; e.g. boot_ngc6791, tuc_grid_nearest_bin02)": dict(base),
    "47tuc, written by repo code (split_ast true; tuc_grid_nearest_bin02_split)": {**base, "split_ast": True},
    "written by an intermediate dev copy (resample_ast only; smooth_test/newmodel)": {**base, "resample_ast": False},
}
for name, have in cases.items():
    with tempfile.TemporaryDirectory() as out:
        os.makedirs(os.path.join(out, "bootstrap"))
        json.dump(have, open(os.path.join(out, "bootstrap", "cell.json"), "w"))
        try:
            B.run_slice(ctx, out, CELL, 2, 0, 1)
            print(f"[ACCEPTED] {name}")
        except SystemExit as e:
            print(f"[REFUSED ] {name}\n           {str(e)[-200:]}")
        except TypeError as e:
            print(f"[TypeError] {name}: {e}")
# what a fresh directory records now
with tempfile.TemporaryDirectory() as out:
    try:
        B.run_slice(ctx, out, CELL, 2, 0, 1)
        print("fresh cell.json:", json.load(open(os.path.join(out, "bootstrap", "cell.json"))))
    except Exception as e:
        print("fresh run failed", e)
