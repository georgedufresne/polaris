#!/usr/bin/env python3
"""Write our Monte Carlo-on 47 Tuc DSEP4 isochrones (../tuc_iso_d4v2_mcfix, as used by ../tuc_final_mc) in the block
format M. Ying's MakeTestCMD.f90 reads, as feh077cmd.<run> (the name his TestCMDPAR.sh passes for 47 Tuc:
global_var.sCMD_vars("47Tuc") -> feh '077').

Magnitudes are exactly those our pipeline uses: read_isochrones + select_age(exact) + BCTable.absolute from
../tuc_final_mc/code with ../tuc_final_mc/setup.json's BC table (MIST ACS/WFC, [Fe/H] -0.75). Absolute F606W and
F606W-F814W; MakeTestCMD and his KS add DM and E(V-I) later.

Ages: his grid, 8-16 Gyr in 0.2 (41 ages; runsCMD/TestCMDPAR take the age in Myr), below each run's age limit
(../tuc_final_mc/tuc_age_limits_d4mc.json). Each block ends at the RGB tip (our isochrones stop there).

MakeTestCMD's reader per age block: one skipped line; header (1x,I3,F9.6,F8.4,F7.0) = npts, mixing length, ALPHAC,
age (Myr) (only npts and age are used; the other two are written as 0); one skipped line; npts list-directed rows
"EEP mass F606W F606W-F814W"; two skipped lines; at most 100 blocks. Rows are written with 6 decimals.

Points without magnitudes (outside the BC table: only the coolest upper-RGB points, e.g. run 10076 at 11 Gyr, log g
-0.42, 2,500 K) are dropped; they must lie more than 20 points above the turn-off. A block is skipped, and listed in
skipped_blocks.csv, if it has fewer than 100 points or does not reach 3 mag above or 4 mag below his subgiant point
(deltaMag 3.0 fit range; e.g. run 10076 at 16 Gyr has 11 points).

Checks per block: npts <= 999; mass strictly increasing; the bluest point (his turn-off) is a main-sequence/turn-off
point (log g > 3.5) and lies before the end; every header parses back with his fixed columns; the written rows read
back to within 5e-7 mag of the pipeline's values.
Output: iso/feh077cmd.<run> here, ages_per_run.json, logs/write_martin_isochrones_tuc.out."""
import json
import os
import sys

import numpy as np

C = "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_final_mc"
sys.path.insert(0, f"{C}/code")
from isochronetoolbox.scmd.model import BCTable, read_isochrones, select_age   # noqa: E402

E = "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_martin_replica_mc"
OUT = f"{E}/iso"
os.makedirs(OUT, exist_ok=True)
setup = json.load(open(f"{C}/setup.json"))
bc = BCTable(setup["bc"])
lim = {int(k): v for k, v in json.load(open(f"{C}/tuc_age_limits_d4mc.json")).items()}
AGES = [a / 1000 for a in range(8000, 16001, 200)]


def fortran_header(s):
    """(1x,I3,F9.6,F8.4,F7.0), blanks ignored."""
    return int(s[1:4].replace(" ", "")), float(s[21:28].replace(" ", ""))


runs = [int(x) for x in open(f"{E}/runs.txt").read().split()]
nblk, nfiles, worst = 0, 0, 0.0
per_run = {}
skipped = []
for r in runs:
    iso = read_isochrones(setup["iso"].format(run=r))
    txt, blocks, ages_written = [], [], []
    for a in AGES:
        if a >= lim[r] - 1e-6:
            continue
        age_yr, rows = select_age(iso, a, exact=True)
        assert abs(age_yr / 1e9 - a) < 1e-6, (r, a, age_yr)
        m = bc.absolute(rows)
        v, vi = m["F606W"], m["F606W"] - m["F814W"]
        mass = rows["initial_mass"].values
        logg = rows["log_g"].values
        ok = np.isfinite(v) & np.isfinite(vi)
        if len(rows) < 100:
            skipped.append((r, a, len(rows), "fewer than 100 points"))
            continue
        ito0 = int(np.nanargmin(np.where(ok, vi, np.nan)))
        assert np.all(np.where(~ok)[0] > ito0 + 20), (r, a, "magnitudes missing below the upper RGB")
        v, vi, mass, logg = v[ok], vi[ok], mass[ok], logg[ok]
        n = len(v)
        assert n <= 999, (r, a, n)
        assert np.all(np.diff(mass) > 0), (r, a, "mass not increasing")
        ito = int(np.argmin(vi))
        assert logg[ito] > 3.5 and ito < n - 5, (r, a, ito, n)
        seg = slice(ito, min(ito + 101, n))
        vsgb = np.interp(vi[ito] + 0.05, vi[seg], v[seg])
        if v.min() > vsgb - 3.0 or v[:ito].max() < vsgb + 4.0:
            skipped.append((r, a, n, f"range: Vsgb {vsgb:.2f}, V {v.min():.2f} to {v[:ito].max():.2f}"))
            continue
        hdr = f"#{n:3d}{0.0:9.6f}{0.0:8.4f}{a * 1000:6.0f}."
        assert len(hdr) == 28
        body = [f"{k + 1} {mass[k]:.6f} {v[k]:.6f} {vi[k]:.6f}" for k in range(n)]
        txt += [f"#NPTS {n}", hdr, "# EEP Mass F606W F606W-F814W", *body, "", ""]
        blocks.append((n, a, v, vi, mass))
        ages_written.append(round(a * 1000))
    assert len(blocks) <= 100
    if not blocks:                     # age limit at or below 8 Gyr (e.g. 10188): no usable age, as in our fit
        per_run[r] = []
        continue
    p = f"{OUT}/feh077cmd.{r}"
    open(p, "w").write("\n".join(txt) + "\n")
    L = open(p).read().splitlines()
    i = 0
    for n, a, v, vi, mass in blocks:
        hn, hage = fortran_header(L[i + 1])
        assert hn == n and abs(hage - a * 1000) < 0.01, (r, L[i + 1])
        d = np.array([[float(x) for x in s.split()] for s in L[i + 3:i + 3 + n]])
        worst = max(worst, np.abs(d[:, 2] - v).max(), np.abs(d[:, 3] - vi).max())
        assert np.abs(d[:, 1] - mass).max() < 5e-7
        assert L[i + 3 + n] == "" and L[i + 4 + n] == ""
        i += 5 + n
        nblk += 1
    assert i == len(L)
    per_run[r] = ages_written
    nfiles += 1
json.dump(per_run, open(f"{E}/ages_per_run.json", "w"))
import pandas as pd                                                   # noqa: E402
pd.DataFrame(skipped, columns=["simulation_id", "age_gyr", "npts", "reason"]).to_csv(f"{E}/skipped_blocks.csv", index=False)
print(f"skipped blocks: {len(skipped)} (skipped_blocks.csv)")
n41 = sum(1 for k in per_run.values() if len(k) == 41)
n0 = sum(1 for k in per_run.values() if not k)
print(f"{nfiles} files, {nblk} age blocks in {OUT}; {n41} runs have all 41 ages, the rest stop at their age limit "
      f"(fewest {min(len(k) for k in per_run.values() if k)}); {n0} runs have no usable age and get no file; every header parses back with MakeTestCMD's columns; largest magnitude "
      f"difference from the pipeline's values after writing {worst:.1e}")
