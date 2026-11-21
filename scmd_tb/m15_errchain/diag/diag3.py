#!/usr/bin/env python3
"""Per-star test: SGB F606W offset vs the local upper-MS colour offset (150 nearest MS stars on the detector).
For a cell shift dF606W = k dcolour, the slope is (k - s_sg) / (1 - s_ms k): reddening k = 2.8,
F606W zero point k = 1, F814W zero point k = 0."""
import json, sys
import numpy as np
from scipy.spatial import cKDTree
E = "/dartfs-hpc/scratch/f007058/scmd_tb/m15_errchain"
sys.path.insert(0, f"{E}/code")
from isochronetoolbox.scmd import Setup, context
S = Setup(**{k: v for k, v in json.load(open(f"{E}/setup.json")).items() if k in Setup.__dataclass_fields__})
C = context(S)
fit = C.ast._args[0]; xc, yc = C.ast.centre_xy
v, c, x, y = fit.v.values, (fit.v - fit.i).values, fit.x.values, fit.y.values
r = np.hypot(x - xc, y - yc)
rid = C.rcut.ridge
m1 = rid.mag > 19.3; o1 = np.argsort(rid.mag[m1])
ms_c = lambda V: np.interp(V, rid.mag[m1][o1], rid.color[m1][o1])
m2 = (rid.mag > 18.2) & (rid.mag < 19.0) & (rid.color > 0.52) & (rid.color < 0.66); o2 = np.argsort(rid.color[m2])
sgb_v = lambda col: np.interp(col, rid.color[m2][o2], rid.mag[m2][o2])
h = 0.01
MS = (v > 19.4) & (v < 20.6) & (r > 300)
SG = (c > 0.53) & (c < 0.64) & (v > 18.2) & (v < 19.0) & (np.abs(v - sgb_v(c)) < 0.2) & (r > 400)
s_ms = np.median((ms_c(v[MS] + h) - ms_c(v[MS] - h)) / (2 * h))
s_sg = np.median((sgb_v(c[SG] + h) - sgb_v(c[SG] - h)) / (2 * h))
ms_res = c[MS] - ms_c(v[MS])
tree = cKDTree(np.column_stack([x[MS], y[MS]]))
for K in (100, 150, 300):
    _, idx = tree.query(np.column_stack([x[SG], y[SG]]), k=K)
    dc = np.median(ms_res[idx], axis=1)
    dv = v[SG] - sgb_v(c[SG])
    # robust line: medians in 8 quantile bins of dc, slope by least squares on the bin medians, bootstrap errors
    def fitline(ii):
        q = np.quantile(dc[ii], np.linspace(0, 1, 9))
        b = np.clip(np.digitize(dc[ii], q[1:-1]), 0, 7)
        X = np.array([np.median(dc[ii][b == j]) for j in range(8)]); Y = np.array([np.median(dv[ii][b == j]) for j in range(8)])
        return np.polyfit(X, Y, 1)[0]
    rs = np.random.RandomState(0)
    sl = fitline(np.arange(len(dc)))
    bs = np.array([fitline(rs.randint(0, len(dc), len(dc))) for _ in range(500)])
    noise = 1.25 * 1.4826 * np.median(np.abs(ms_res - np.median(ms_res))) / np.sqrt(K)
    att = np.var(dc) / (np.var(dc) + 0) if False else (np.var(dc) - noise ** 2) / np.var(dc)
    print(f"K {K}: {SG.sum()} SGB stars, local colour offset rms {np.std(dc):.4f} (noise ~{noise:.4f}); "
          f"slope {sl:.2f} +/- {bs.std():.2f}; corrected for noise in dc: {sl / att:.2f}")
print(f"ridge slopes: MS dcol/dV {s_ms:.4f}, SGB dV/dcol {s_sg:.2f}; predicted slopes:")
for name, k in (("reddening k=2.8", 2.8), ("F606W zero point k=1", 1.0), ("F814W zero point k=0", 0.0)):
    print(f"  {name:22s}: {(k - s_sg) / (1 - s_ms * k):.2f}")
