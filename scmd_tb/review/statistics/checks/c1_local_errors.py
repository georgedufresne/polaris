"""Check 1: real 47 Tuc ASTs -- local error distributions (k=30, crowding metric) and copy vs smooth draws."""
import sys, numpy as np, pandas as pd
sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/review/statistics/dev")
from isochronetoolbox.scmd.ast import NearestErrors, centre
from scipy import stats
F32 = np.float32
A = pd.read_csv("/dartfs-hpc/scratch/f007058/Files_from_Martin/47Tuc_data/47Tucartstars.dat", dtype=F32)
zpv, zpi = 30.5742, 29.6777
split = (16.48 - zpv, 15.55 - zpi)
print("n AST", len(A))
ne = NearestErrors(A, k=30, smooth=True, metric="crowding", split=split)
ne100 = NearestErrors(A, k=100, smooth=True)            # linear metric, k=100 (smooth_resample config)
x, y = A.InputX.values, A.InputY.values
xc, yc = centre(x, y); r = np.hypot(x - xc, y - yc)
vin = A.InputF606W.values + zpv
print("AST centre", xc, yc, " recovered frac", ne.rec.mean())
# completeness and error scale by radius
for lo, hi in [(0, 100), (100, 200), (200, 400), (400, 800), (800, 1500), (1500, 4000)]:
    s = (r >= lo) & (r < hi) & (vin > 15.0) & (vin < 19.1)
    rec = ne.rec[s]
    ev = ne.ev[s][rec]; ei = ne.ei[s][rec]
    mad = 1.4826 * np.median(np.abs(ev - np.median(ev)))
    print(f"r {lo:4d}-{hi:4d}: n {s.sum():6d} compl {rec.mean():.3f}  dV MAD-sigma {mad:.4f} SD {ev.std():.4f} "
          f"SD/MADsig {ev.std()/mad:.2f} kurt(excess) {stats.kurtosis(ev):6.1f} skew {stats.skew(ev):5.2f} "
          f"mean-median {ev.mean()-np.median(ev):+.4f} corr(dV,dI) {np.corrcoef(ev, ei)[0,1]:.2f}")
# neighbourhood statistics: recovered count n, h, c
nb = ne.nbr; n = ne.rec[nb].sum(1)
h = 0.5 * np.maximum(n, 1) ** (-1/6); c = 1/np.sqrt(1 + h**2)
print("recovered neighbours n: pct 1,10,50,90:", np.percentile(n, [1, 10, 50, 90]), " h median", np.median(h), " c min/median", c.min(), np.median(c))
print("MAD floor hits (sv<=1e-4):", (ne._L[:,0] <= 0.5 * 1e-4 * 1.001).sum())
# per-neighbourhood: within-neighbourhood SD / robust sigma and mean-median of recovered values
ok = ne.rec[nb]
evn = np.where(ok, ne.ev[nb], np.nan)
sd = np.nanstd(evn, 1); med = np.nanmedian(evn, 1); mn = np.nanmean(evn, 1)
madsig = 1.4826 * np.nanmedian(np.abs(evn - med[:, None]), 1)
q = (n >= 10)
ratio = (sd / madsig)[q]
print("per-neighbourhood SD/MAD-sigma pct 10,50,90,99:", np.round(np.nanpercentile(ratio, [10, 50, 90, 99]), 2))
# implied variance change of smooth draws: Var(Y)/Var(X) = (Var + h^2 S_r)/((1+h^2) Var)
hh = h[q]
vr = (sd[q]**2 + hh**2 * madsig[q]**2) / ((1 + hh**2) * sd[q]**2)
print("smooth: within-neighbourhood Var ratio pct 1,50,99:", np.round(np.nanpercentile(vr, [1, 50, 99]), 3))
mshift = ((1 - c[q]) * (mn - med)[q])
print("smooth: mean shift toward median |(1-c)(mean-median)| pct 50,90,99 (mag):", np.round(np.nanpercentile(np.abs(mshift), [50, 90, 99]), 5))
# aggregated draws in radial zones: copy vs smooth, same random numbers
rs0 = np.random.RandomState(1)
for lo, hi in [(0, 150), (150, 400), (400, 1000), (1000, 3000)]:
    m = 200000
    rr = np.sqrt(rs0.uniform(lo**2, hi**2, m))
    vv = rs0.uniform(16.0, 18.9, m)
    # place on the sequence: i from the ASTs' own (V_in, I_in) relation
    sel = (vin > 15.5) & (vin < 19.5)
    pf = np.polyfit(vin[sel], A.InputF814W.values[sel] + zpi, 3)
    ii = np.polyval(pf, vv)
    args = (rr, vv - zpv, ii - zpi)
    neC = ne; ok_s, ev_s, ei_s = neC.draw(np.random.RandomState(5), *args)
    neC.smooth = False; ok_c, ev_c, ei_c = neC.draw(np.random.RandomState(5), *args); neC.smooth = True
    a, b = ev_c[ok_c], ev_s[ok_s]
    mad = lambda z: 1.4826 * np.median(np.abs(z - np.median(z)))
    col_c, col_s = (ev_c - ei_c)[ok_c], (ev_s - ei_s)[ok_s]
    tail = lambda z, t: np.mean(np.abs(z - np.median(z)) > t)
    print(f"zone r {lo}-{hi}: dV MAD copy {mad(a):.4f} smooth {mad(b):.4f} ({mad(b)/mad(a)-1:+.3f}); SD {a.std():.4f} {b.std():.4f} ({b.std()/a.std()-1:+.3f}); "
          f"colour MAD {mad(col_c):.4f} {mad(col_s):.4f} ({mad(col_s)/mad(col_c)-1:+.3f}) SD {col_c.std():.4f} {col_s.std():.4f}; "
          f"|dV-med|>0.1: {tail(a,0.1):.4f} {tail(b,0.1):.4f}; >0.3: {tail(a,0.3):.5f} {tail(b,0.3):.5f}; mean {a.mean():+.5f} {b.mean():+.5f}")
