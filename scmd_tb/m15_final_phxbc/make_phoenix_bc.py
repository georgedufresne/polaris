#!/usr/bin/env python3
"""Write Dartmouth PHOENIX ACS/WFC bolometric corrections in the format of the pipeline's BCTable (npz: points (Teff,
log g), bc per band with M = M_bol - bc, feh, bands), for the one-change test against production's MIST table.

Source: /dartfs-hpc/rc/lab/C/ChaboyerB/Catherine/dsep3/iso/phx/Zm2d0ap4.HST_ACSWF and Zm2d5ap4.HST_ACSWF (columns Teff,
log g, [Fe/H], BCv, V-F435W, V-F475W, V-F555W, V-F606W, V-F625W, V-F775W, V-F814W, V-F850LP; BCv = M_V - M_bol here, the
Sun's entry being +0.082). [Fe/H] -2.25 is the mean of the -2.0 and -2.5 tables on the (Teff, log g) points the two
share, matching production's single MIST table at -2.25. [alpha/Fe] = +0.4 is the most common value in the runs
(0.2 vs 0.4 changes the colour by <= 0.004). bc_X = -BCv + (V - X), so that M_X = M_bol - bc_X = M_bol + BCv - (V - X).

Check: on the isochrones of runs 10357, 10392 and 10089 (8-16 Gyr), BCTable with this file is compared with
../m15_isochrone_bc_test's PHOENIX conversion at [Fe/H] -2.25, [alpha/Fe] 0.4 (the maximum difference is printed).
Output: inputs/bc_phoenix_acs_wfc_fehm225_ap4.npz, logs/make_phoenix_bc.out."""
import sys

import numpy as np
import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
O = f"{B}/m15_final_phxbc"
P = "/dartfs-hpc/rc/lab/C/ChaboyerB/Catherine/dsep3/iso/phx"
cols = ["teff", "logg", "feh", "bcv", "v435", "v475", "v555", "v606", "v625", "v775", "v814", "v850"]
a = pd.read_csv(f"{P}/Zm2d0ap4.HST_ACSWF", sep=r"\s+", comment="#", header=None, names=cols)
b = pd.read_csv(f"{P}/Zm2d5ap4.HST_ACSWF", sep=r"\s+", comment="#", header=None, names=cols)
m = a.merge(b, on=["teff", "logg"], suffixes=("_a", "_b"))
print(f"-2.0 table {len(a)} points, -2.5 table {len(b)}, shared {len(m)}")
bc606 = 0.5 * ((-m.bcv_a + m.v606_a) + (-m.bcv_b + m.v606_b))
bc814 = 0.5 * ((-m.bcv_a + m.v814_a) + (-m.bcv_b + m.v814_b))
out = f"{O}/inputs/bc_phoenix_acs_wfc_fehm225_ap4.npz"
np.savez(out, points=m[["teff", "logg"]].values.astype(float), bc=np.column_stack([bc606, bc814]),
         bands=np.array(["F606W", "F814W"]), feh=np.float64(-2.25), av=np.float64(0.0), rv=np.float64(3.1))
print(f"wrote {out}")

sys.path.insert(0, f"{B}/m15_final/code")
from isochronetoolbox.scmd.model import BCTable, read_isochrones    # noqa: E402
sys.argv = ["x"]
src = open(f"{B}/m15_isochrone_bc_test/bc_test.py").read()
ns = {}
exec(src[src.index("class Phoenix"):src.index("phx = Phoenix()")], {"np": np, "pd": pd, "LinearNDInterpolator":
     __import__("scipy.interpolate", fromlist=["LinearNDInterpolator"]).LinearNDInterpolator, "PHX": P}, ns)
phx = ns["Phoenix"]()
t = BCTable(out)
for run in (10357, 10392, 10089):
    iso = read_isochrones(f"/dartfs-hpc/scratch/f007058/M15/iso_dsep4/{run}/isos/isochrones_example.txt")
    iso = iso[(iso.isochrone_age_yr >= 8e9) & (iso.isochrone_age_yr <= 16e9)]
    x, y = t.absolute(iso), phx.absolute(iso, -2.25, 0.4)
    win = np.isfinite(y["F606W"]) & (y["F606W"] > 0.5) & (y["F606W"] < 7.5)       # absolute range of the fit and calibration
    d6, dc = np.abs(x["F606W"] - y["F606W"])[win], np.abs((x["F606W"] - x["F814W"]) - (y["F606W"] - y["F814W"]))[win]
    print(f"run {run}: {win.sum()} points with M_F606W 0.5-7.5; max |dF606W| {np.nanmax(d6):.2e}, max |d colour| "
          f"{np.nanmax(dc):.2e}, 99th pct |d colour| {np.nanpercentile(dc, 99):.2e}")
