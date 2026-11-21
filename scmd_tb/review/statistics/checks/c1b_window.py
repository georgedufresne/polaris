"""Check 1b: restrict to the fit window; split by neighbourhood size; analytic check on a Laplace/heavy-tail model."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/review/statistics/dev")
from isochronetoolbox.scmd.ast import NearestErrors, centre
F32 = np.float32
A = pd.read_csv("/dartfs-hpc/scratch/f007058/Files_from_Martin/47Tuc_data/47Tucartstars.dat", dtype=F32)
zpv, zpi = 30.5742, 29.6777
split = (16.48 - zpv, 15.55 - zpi)
mad = lambda z: 1.4826 * np.median(np.abs(z - np.median(z)))
for k, metric, sp in [(30, "crowding", split), (100, "linear", None)]:
    ne = NearestErrors(A, k=k, smooth=True, metric=metric, split=sp)
    vin = A.InputF606W.values + zpv
    win = (vin > 15.02) & (vin < 19.02)
    n = ne.rec[ne.nbr].sum(1)
    print(f"k={k} {metric}: in-window recovered neighbours pct 1,5,10,50: {np.percentile(n[win], [1, 5, 10, 50])}; "
          f"floor sv=1e-4 in window: {(ne._L[win,0] <= ne.bw * max(n[win].max(),1)**(-1/6) * 1.0001e-4 + 1e-12).sum()}")
    # draw at the ASTs' own input positions in the window (realistic coverage), paired copy vs smooth
    x, y = A.InputX.values, A.InputY.values; xc, yc = centre(x, y); r = np.hypot(x - xc, y - yc)
    idx = np.flatnonzero(win)
    rs = np.random.RandomState(0); idx = rs.choice(idx, 300000)
    args = (r[idx], A.InputF606W.values[idx].astype(float), A.InputF814W.values[idx].astype(float))
    ok_s, ev_s, ei_s = ne.draw(np.random.RandomState(9), *args)
    ne.smooth = False; ok_c, ev_c, ei_c = ne.draw(np.random.RandomState(9), *args); ne.smooth = True
    _, near = ne.tree.query(ne._coords(args[2-2+1], args[2], args[0]))
    nn = n[near]
    for lo, hi in [(0, 200), (200, 600), (600, 1500), (1500, 5000)]:
        for nlo in (0, 20):
            s = ok_c & (r[idx] >= lo) & (r[idx] < hi) & (nn >= nlo)
            a, b = ev_c[s], ev_s[s]; ca, cb = (ev_c - ei_c)[s], (ev_s - ei_s)[s]
            t3 = lambda z: np.mean(np.abs(z - np.median(z)) > 0.3)
            print(f"   r {lo:4d}-{hi:4d} n>={nlo:2d}: N {s.sum():6d}  dV MAD {mad(b)/mad(a)-1:+.3f}  SD {b.std()/a.std()-1:+.3f}  "
                  f"colour MAD {mad(cb)/mad(ca)-1:+.3f} SD {cb.std()/ca.std()-1:+.3f}  tail>0.3 {t3(a):.4f}->{t3(b):.4f}")
# analytic: one neighbourhood = n draws from a cusp + heavy tail; smooth exactly as _kernel/draw do
rs = np.random.RandomState(3)
def smooth_once(e, bw=0.5):
    m = np.median(e); s = 1.4826 * np.median(np.abs(e - m)); h = bw * len(e) ** (-1/6); c = 1/np.sqrt(1 + h*h)
    return m, s, h, c
for name, gen in [("Gaussian", lambda n: rs.normal(0, 1, n)), ("Laplace", lambda n: rs.laplace(0, 1, n)),
                  ("Laplace+5% exp tail", lambda n: rs.laplace(0, 1, n) + np.where(rs.random_sample(n) < 0.05, rs.exponential(8, n), 0))]:
    for bw in (1.0, 0.5, 0.25):
        rc, rsd, rt, ra = [], [], [], []
        for rep in range(4000):
            e = gen(30)
            m, s, h, c = smooth_once(e, bw)
            pick = e[rs.randint(0, 30, 200)]
            y = m + c * (pick - m + h * s * rs.normal(0, 1, 200))
            rc.append(np.c_[pick, y])
        P = np.vstack(rc)
        a, b = P[:, 0], P[:, 1]
        print(f"{name:20s} bw {bw}: MAD {mad(b)/mad(a)-1:+.3f}  SD {b.std()/a.std()-1:+.3f}  P(|x|>5) {np.mean(np.abs(a)>5):.4f}->{np.mean(np.abs(b)>5):.4f}")
