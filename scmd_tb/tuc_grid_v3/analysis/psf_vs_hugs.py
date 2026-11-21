#!/usr/bin/env python3
"""Star-by-star scatter between the Treasury and HUGS reductions of the same F606W/F814W
images, per F606W bin, against the PSF scatter injected in v2."""
import json
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
B = "/dartfs-hpc/scratch/f007058/scmd_tb"
acs = pd.read_csv(json.load(open(f"{B}/tuc_grid_v2/setup.json"))["obs_raw"], sep=r"\s+")
hug = pd.read_csv(f"{B}/hugs/47Tuc_hugs_meth2.dat", sep=r"\s+")
tab = np.array(json.load(open(f"{B}/tuc_grid_v2/setup.json"))["extra_sigma_table"])
d, j = cKDTree(hug[["x", "y"]].values).query(acs[["x", "y"]].values)
ok = (d < 0.5) & (np.abs(acs.Vvega.values - hug.Vvega.values[j]) < 0.1)
v = acs.Vvega.values
dv = acs.Vvega.values - hug.Vvega.values[j]
di = acs.Ivega.values - hug.Ivega.values[j]
dc = acs.VIvega.values - hug.VIvega.values[j]
# keep stars on the sequence: |colour - local median| small in the Treasury catalogue
mad = lambda x: 1.4826 * np.median(np.abs(x - np.median(x)))
print(f"{'F606W bin':>12} {'N':>6} {'dF606W':>8} {'dF814W':>8} {'dcolour':>8} {'r(dV,dI)':>9} | {'PSF sigma/band':>14} {'PSF colour':>10}")
for lo, hi in [(15.5, 16.4), (17.6, 18.0), (18.0, 18.4), (18.4, 18.8), (18.8, 19.0)]:
    m = ok & (v > lo) & (v < hi)
    a, b, c = dv[m] - np.median(dv[m]), di[m] - np.median(di[m]), dc[m] - np.median(dc[m])
    core = (np.abs(a) < 5 * mad(a)) & (np.abs(b) < 5 * mad(b))
    r = np.corrcoef(a[core], b[core])[0, 1]
    s = np.interp((lo + hi) / 2, tab[:, 0], tab[:, 1])
    print(f"{lo:5.1f}-{hi:4.1f}  {m.sum():6d} {mad(a)*1e3:8.1f} {mad(b)*1e3:8.1f} {mad(c)*1e3:8.1f} {r:9.2f} | "
          f"{s*1e3:14.1f} {np.sqrt(2)*s*1e3:10.1f}")
print("(mmag; robust sigma; per catalogue, if the two reductions' own errors are equal and independent, divide by sqrt 2)")
