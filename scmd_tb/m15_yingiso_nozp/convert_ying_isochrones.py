#!/usr/bin/env python3
"""Write Ying et al. (2025)'s published M15 isochrones (Zenodo 10.5281/zenodo.13998637, M15_MC_iso.h5) for our 500 runs
in the DSEP block format the pipeline already reads (model.read_dsep, iso_format "dsep"): per age a "#NPTS" line, a
"#<npts> <mixlen> <ALPHAC> <age Myr>." line, a header, then rows EEP, Mass, F606W, F606W-F814W (absolute, his values,
written with full precision). His file has 41 ages, 8-16 Gyr in 0.2 Gyr steps.

Also writes age_limits_ying.json: per run, the first age whose isochrone does not span absolute F606W < 1.0 to > 6.0
(the fit window is F606W 16.64-20.64, i.e. absolute 1.0-5.3 for mu 15.3-15.6), else 16.2.

Validation: every file is read back with this experiment's read_dsep and compared with the HDF5 values (exact).
Output: /dartfs-hpc/scratch/f007058/M15/iso_ying2025/<run>/ying_cmd.<run>, age_limits_ying.json there,
logs/convert_ying_isochrones.out."""
import json
import os
import sys

import h5py
import numpy as np
import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
H = f"{B}/m15_yingiso_nozp"
OUT = "/dartfs-hpc/scratch/f007058/M15/iso_ying2025"
sys.path.insert(0, f"{H}/code")
from isochronetoolbox.scmd.model import read_dsep   # noqa: E402

YF = h5py.File("/dartfs-hpc/scratch/f007058/M15/zenodo_ying2025/M15_MC_iso.h5", "r")
runs = [int(x) for x in open(f"{H}/runs.txt").read().split()]


def frame(g):
    return pd.DataFrame(g["block0_values"][()], columns=[x.decode() for x in g["block0_items"][()]])


lim, bad, nmax = {}, 0, 0.0
for run in runs:
    var = frame(YF[f"{run}/var"]).iloc[0]
    cmd = YF[f"{run}/isochrones/cmd"]
    ages = sorted(int(a) for a in cmd)
    blocks, ref, lim[run] = [], [], 16.2
    for a in ages:
        t = frame(cmd[str(a)]).sort_values("EEP")
        if lim[run] > 16.1 and not (t.F606W.min() < 1.0 and t.F606W.max() > 6.0):
            lim[run] = a / 1000
        n = len(t)
        body = "\n".join(f"{int(e)} {m!r} {v!r} {c!r}" for e, m, v, c in zip(t.EEP, t.Mass, t.F606W, t.F814W))
        blocks.append(f"#NPTS {n}\n#{n} {var.CMIXLA!r} {var.ALPHAC!r} {a}.\n# EEP Mass F606W F606W-F814W\n{body}")
        ref.append(pd.DataFrame({"isochrone_age_yr": a * 1e6, "initial_mass": t.Mass.values,
                                 "F606W": t.F606W.values, "F814W": t.F606W.values - t.F814W.values}))
    os.makedirs(f"{OUT}/{run}", exist_ok=True)
    path = f"{OUT}/{run}/ying_cmd.{run}"
    open(path, "w").write("\n".join(blocks) + "\n")
    back = read_dsep(path)
    ref = pd.concat(ref, ignore_index=True)
    ok = len(back) == len(ref) and all(np.array_equal(back[c].values, ref[c].values) for c in ref.columns)
    nmax = max(nmax, float(np.max(np.abs(back[ref.columns].values - ref.values))) if len(back) == len(ref) else np.inf)
    bad += not ok
json.dump({str(k): v for k, v in lim.items()}, open(f"{OUT}/age_limits_ying.json", "w"), indent=0)
v = np.array(list(lim.values()))
print(f"{len(runs)} runs written to {OUT}; read back with read_dsep: {len(runs) - bad} exact, {bad} not "
      f"(max abs difference {nmax:g}); ages per run: {len(ages)} ({ages[0] / 1000}-{ages[-1] / 1000} Gyr)")
print(f"age limits: {int((v > 16.1).sum())} runs valid to 16 Gyr; the others first fail at "
      f"{sorted(set(v[v < 16.1]))[:10]}")
