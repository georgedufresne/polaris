#!/usr/bin/env python3
"""For the 1000-run repeat: write Ying's isochrones for runs 10500-10999 (beyond our first 500) in MakeTestCMD's block
format, straight from his HDF5 (M15_MC_iso.h5, Zenodo 10.5281/zenodo.13998637), as write_martin_isochrones.py did
from the ying_cmd copies; and compute calibration-star weights for all 1000 runs with the code of
../m15_ying_isochrones/ying_calib.py (HIP 46120 and HIP 106924, normalised point-to-isochrone distance, minimum over
8-16 Gyr, p = exp(-chi2/2)).

Rows are written as in ../m15_yingiso_nozp/convert_ying_isochrones.py: EEP, Mass, F606W, F606W-F814W at full precision
(repr). Validation: headers parse back with MakeTestCMD's fixed columns, rows read back equal to the HDF5 values, and
for runs 10000-10499 the rows equal the files already written (same writer, same source) and the calibration chi2
equals ../m15_ying_isochrones/ying_calib_chi2.csv.
Output: M15/iso_ying2025_martinfmt/feh227cmd.<run>, runs_1000.txt, ying_calib_chi2_1000.csv,
logs/write_more_isochrones.out."""
import h5py
import numpy as np
import pandas as pd

E = "/dartfs-hpc/scratch/f007058/scmd_tb/m15_martin_replica"
OUT = "/dartfs-hpc/scratch/f007058/M15/iso_ying2025_martinfmt"
YF = h5py.File("/dartfs-hpc/scratch/f007058/M15/zenodo_ying2025/M15_MC_iso.h5", "r")
STARS = {"HIP 46120": (5.7867, 0.0026, 0.566, 0.002), "HIP 106924": (6.0406, 0.0037, 0.601, 0.005)}
RUNS = list(range(10000, 11000))


def frame(g):
    return pd.DataFrame(g["block0_values"][()], columns=[x.decode() for x in g["block0_items"][()]])


def fortran_header(s):
    return (int(s[1:4].replace(" ", "")), float(s[4:13].replace(" ", "")), float(s[13:21].replace(" ", "")),
            float(s[21:28].replace(" ", "")))


def dist2(x, y, x0, y0):
    ax, ay, bx, by = x[:-1], y[:-1], x[1:], y[1:]
    dx, dy = bx - ax, by - ay
    t = np.clip(((x0 - ax) * dx + (y0 - ay) * dy) / np.maximum(dx * dx + dy * dy, 1e-30), 0, 1)
    return np.min((ax + t * dx - x0) ** 2 + (ay + t * dy - y0) ** 2)


cal, same_rows = [], 0
for run in RUNS:
    var = frame(YF[f"{run}/var"]).iloc[0]
    cmd = YF[f"{run}/isochrones/cmd"]
    txt, best = [], None
    for a in sorted(int(x) for x in cmd):
        t = frame(cmd[str(a)]).sort_values("EEP")
        n = len(t)
        assert n <= 999
        rows = [f"{int(e)} {m!r} {v!r} {c!r}" for e, m, v, c in zip(t.EEP, t.Mass, t.F606W, t.F814W)]
        hdr = f"#{n:3d}{var.CMIXLA:9.6f}{var.ALPHAC:8.4f}{a:6.0f}."
        assert len(hdr) == 28
        hn, hml, _, hage = fortran_header(hdr)
        assert hn == n and abs(hage - a) < 0.01 and abs(hml - var.CMIXLA) < 1e-6
        back = np.array([[float(x) for x in r.split()] for r in rows])
        assert np.array_equal(back, t[["EEP", "Mass", "F606W", "F814W"]].values)
        txt += [f"#NPTS {n}", hdr, "# EEP Mass F606W F606W-F814W", *rows, "", ""]
        if 8000 <= a <= 16000:                                 # calibration, as ying_calib.py
            v, c = t.F606W.values, t.F814W.values
            chi2 = sum(dist2(c / sC, v / sM, C / sC, M / sM) for M, sM, C, sC in STARS.values())
            if best is None or chi2 < best[1]:
                best = (a / 1000, chi2)
    p = f"{OUT}/feh227cmd.{run}"
    if run < 10500:
        same_rows += open(p).read() == "\n".join(txt) + "\n"
    else:
        open(p, "w").write("\n".join(txt) + "\n")
    cal.append(dict(simulation_id=run, chi2_min=best[1], age_at_min=best[0], p=np.exp(-best[1] / 2)))
c = pd.DataFrame(cal).set_index("simulation_id")
c.to_csv(f"{E}/ying_calib_chi2_1000.csv")
open(f"{E}/runs_1000.txt", "w").write("\n".join(str(r) for r in RUNS) + "\n")
old = pd.read_csv("/dartfs-hpc/scratch/f007058/scmd_tb/m15_ying_isochrones/ying_calib_chi2.csv").set_index("simulation_id")
dmax = np.max(np.abs(c.loc[old.index, "chi2_min"].values - old.chi2_min.values))
print(f"runs 10500-10999: 500 files written to {OUT} (headers parse back, rows equal the HDF5 values)")
print(f"runs 10000-10499: {same_rows} of 500 existing files identical to this writer's output")
print(f"calibration: 1000 runs; max |chi2 - ying_calib_chi2.csv| over the first 500: {dmax:.2e}; chi2 < 2: "
      f"first 500 {(c.chi2_min.iloc[:500] < 2).mean():.0%}, new 500 {(c.chi2_min.iloc[500:] < 2).mean():.0%}")
