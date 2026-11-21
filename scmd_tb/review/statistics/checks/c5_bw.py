"""Check 5: bandwidth sweep on real ASTs (k=30 crowding): core width, SD, far tail vs copy."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/review/statistics/dev")
from isochronetoolbox.scmd.ast import NearestErrors, centre
F32 = np.float32
A = pd.read_csv("/dartfs-hpc/scratch/f007058/Files_from_Martin/47Tuc_data/47Tucartstars.dat", dtype=F32)
zpv, zpi = 30.5742, 29.6777
ne = NearestErrors(A, k=30, smooth=True, metric="crowding", split=(16.48 - zpv, 15.55 - zpi))
x, y = A.InputX.values, A.InputY.values; xc, yc = centre(x, y); r = np.hypot(x - xc, y - yc)
vin = A.InputF606W.values + zpv; win = (vin > 15.02) & (vin < 19.02)
idx = np.random.RandomState(0).choice(np.flatnonzero(win), 300000)
args = (r[idx], A.InputF606W.values[idx].astype(float), A.InputF814W.values[idx].astype(float))
mad = lambda z: 1.4826 * np.median(np.abs(z - np.median(z)))
t3 = lambda z: np.mean(np.abs(z - np.median(z)) > 0.3)
ne.smooth = False; ok_c, ev_c, ei_c = ne.draw(np.random.RandomState(9), *args); ne.smooth = True
zones = [(0, 200), (200, 600), (600, 1500)]
for bw, shrink in [(0.5, True), (0.25, True), (0.15, True), (0.25, False), (0.15, False)]:
    ne.bw = bw; ne._kernel()
    if not shrink:
        ne._c = np.ones_like(ne._c)
    ok_s, ev_s, ei_s = ne.draw(np.random.RandomState(9), *args)
    out = []
    for lo, hi in zones:
        s = ok_c & (r[idx] >= lo) & (r[idx] < hi)
        a, b = ev_c[s], ev_s[s]; ca, cb = (ev_c - ei_c)[s], (ev_s - ei_s)[s]
        out.append(f"r{lo}-{hi}: dV MAD {mad(b)/mad(a)-1:+.3f} col MAD {mad(cb)/mad(ca)-1:+.3f} SD {b.std()/a.std()-1:+.3f} tail {t3(b)/t3(a)-1:+.3f}")
    print(f"bw {bw:4.2f} shrink {str(shrink):5s} | " + " | ".join(out))
    print(f"      unique dV among recovered draws: {len(np.unique(ev_s[ok_s]))} of {ok_s.sum()}")
