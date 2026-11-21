"""Check 7: sensitivity of the 1/c weights to evaluating c at observed rather than true magnitudes."""
import sys, json, numpy as np, pandas as pd
sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/review/statistics/dev")
from isochronetoolbox.scmd import Setup, context
from isochronetoolbox.scmd.ast import centre
F32 = np.float32
S = Setup(**{k: v for k, v in json.load(open("/dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/newmodel/setup.json")).items()
             if k in Setup.__dataclass_fields__})
ctx = context(S); ne = ctx.ast.nearest; zv, zi = S.zp_v, S.zp_i
raw = pd.read_csv(S.obs_raw, sep=r"\s+")
d = raw[ctx.rcut.select(raw.VIvega.values, raw.Vvega.values, ctx.rcut.ridge)]
v, i = d.Vvega.values.astype(float), d.Ivega.values.astype(float)
xc, yc = centre(d.x.values.astype(F32), d.y.values.astype(F32))
R = np.hypot(d.x.values.astype(F32) - F32(xc), d.y.values.astype(F32) - F32(yc)).astype(float)
c0 = ne.completeness(v - zv, i - zi, R)
for dm in (0.03, 0.3):
    c1 = ne.completeness(v + dm - zv, i + dm - zi, R)          # the true star fainter by dm (a blend / an error)
    for lo, hi in [(0, 400), (400, 5000)]:
        s = (R >= lo) & (R < hi)
        print(f"shift +{dm} mag, r {lo}-{hi}: mean |1/c(true) - 1/c(obs)| {np.mean(np.abs(1/c1[s]-1/c0[s])):.4f}, mean 1/c {np.mean(1/c0[s]):.3f}")
for lo, hi in [(0, 400)]:
    for b0 in np.arange(15.0, 19.0, 1.0):
        s = (R >= lo) & (R < hi) & (v >= b0) & (v < b0 + 1)
        print(f"  r {lo}-{hi} V {b0:.0f}-{b0+1:.0f}: mean c {c0[s].mean():.3f} (n {s.sum()})")
