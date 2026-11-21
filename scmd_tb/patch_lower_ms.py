#!/usr/bin/env python3
"""Fill the main-sequence holes of a pyeep isochrone set with the tracks themselves.

pyeep (v1.0 EEP rules) leaves out masses whose track never reaches the TAMS:
on the 47 Tuc tracks the 0.60-0.75 Msun models are still burning hydrogen at
18 Gyr, so no isochrone point is placed for them past the IAMS and a 12 Gyr
isochrone jumps from 0.45 to 0.80 Msun (V ~ 22 -> 18 at 47 Tuc's distance).

Main-sequence stars at these ages are unevolved, so the isochrone there is
just each track's state at that age: for every age block, each track that
(1) reaches the age, (2) is still on the main sequence there (central H > 1e-3)
and (3) falls in a mass gap of the isochrone wider than ``--min-gap`` gets one
row, interpolated linearly in age along the track.  Existing isochrone points
are never replaced.  The result has the same columns; added rows have EEP = 0.

    patch_lower_ms.py ISOCHRONE_FILE MESA_TRACK_DIR OUT_FILE [--min-gap 0.04]
"""
import argparse
import glob
import os

import numpy as np
import pandas as pd

COLS = ["EEP", "isochrone_age_yr", "initial_mass", "star_mass", "star_mdot", "he_core_mass",
        "c_core_mass", "log_L", "log_LH", "log_LHe", "log_Teff", "log_R", "log_g", "surface_h1",
        "surface_he3", "surface_he4", "surface_c12", "surface_o16", "log_center_T",
        "log_center_Rho", "center_gamma", "center_h1", "center_he4", "center_c12"]


def read_track(path):
    lines = open(path).read().split("\n")
    h = next(i for i, l in enumerate(lines) if "star_age" in l.split())
    cols = lines[h].split()
    rows = [l.split() for l in lines[h + 1:] if len(l.split()) == len(cols)]
    return pd.DataFrame(np.array(rows, dtype=float), columns=cols)


def read_blocks(path):
    """Header lines and the data rows of each age block (all 24 columns)."""
    lines = open(path).read().split("\n")
    head = lines[:11]
    rows = [l.split() for l in lines[11:] if l.strip() and not l.lstrip().startswith("#")]
    df = pd.DataFrame(np.array(rows, dtype=float), columns=COLS)
    return head, df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("iso")
    ap.add_argument("tracks")
    ap.add_argument("out")
    ap.add_argument("--min-gap", type=float, default=0.04)
    a = ap.parse_args()
    head, iso = read_blocks(a.iso)
    tracks = []
    for f in sorted(glob.glob(os.path.join(a.tracks, "m[0-9]*"))):
        t = read_track(f)
        tracks.append((float(t["star_mass"].iloc[0]), t))
    out, added = [], 0
    for age, blk in iso.groupby("isochrone_age_yr", sort=True):
        blk = blk.sort_values("initial_mass")
        m = blk["initial_mass"].values
        new = []
        for M, t in tracks:
            age_col = t["star_age"].values
            if age_col[-1] < age or age_col[0] > age:
                continue
            k = np.searchsorted(m, M)
            lo = m[k - 1] if k > 0 else -np.inf
            hi = m[k] if k < len(m) else np.inf
            if not (hi - lo > a.min_gap and min(M - lo, hi - M) > 0.01):
                continue
            row = {c: (np.interp(age, age_col, t[c].values) if c in t.columns else 0.0) for c in COLS}
            if row["center_h1"] <= 1e-3:
                continue
            row.update(EEP=0.0, isochrone_age_yr=age, initial_mass=M)
            new.append(row)
        if new:
            blk = pd.concat([blk, pd.DataFrame(new)[COLS]]).sort_values("initial_mass")
            added += len(new)
        out.append(blk)
    res = pd.concat(out)
    with open(a.out, "w") as fh:
        fh.write("\n".join(head) + "\n")
        for _, blk in res.groupby("isochrone_age_yr", sort=True):
            np.savetxt(fh, blk[COLS].values, fmt="%.16E")
    print(f"{os.path.basename(a.iso)}: {added} track points added over {res.isochrone_age_yr.nunique()} ages")


if __name__ == "__main__":
    main()
