"""Per neighbourhood: what the smooth draw does to the core (MAD) and to the tail (|dev| > 5 MAD)
relative to the neighbourhood's own copies, averaged over the neighbourhoods the sCMD uses."""
import sys; sys.path.insert(0, ".")
from common import *
S = load("dev")
s = setup_fresh(S); ctx = S.context(s); ne = ctx.ast.nearest
rs = np.random.RandomState(0)
# neighbourhoods with many recovered members, inside the magnitude range used
n = ne.rec[ne.nbr].sum(axis=1)
use = np.flatnonzero(n >= 25)[::50]
mad = lambda x: 1.4826 * np.median(np.abs(x - np.median(x)))
r_mad, r_tail, r_sd = [], [], []
for a in use:
    nb = ne.nbr[a][ne.rec[ne.nbr[a]]]
    ev = ne.ev[nb]
    m, L, c = ne._m[a, 0], ne._L[a, 0], ne._c[a]
    M = 4000
    pick = ev[rs.randint(0, len(ev), M)]
    sm = m + c * (pick - m + L * rs.standard_normal(M))
    r_mad.append(mad(sm) / mad(pick)); r_sd.append(sm.std() / pick.std())
    t = 5 * mad(ev)
    r_tail.append((np.abs(sm - m) > t).mean() / max((np.abs(pick - m) > t).mean(), 1e-9) if (np.abs(ev - m) > t).any() else np.nan)
print(f"{len(use)} neighbourhoods: smooth/copy MAD ratio median {np.median(r_mad):.3f} (IQR {np.percentile(r_mad,25):.3f}-{np.percentile(r_mad,75):.3f}); "
      f"SD ratio median {np.median(r_sd):.3f}; P(|dev|>5 MAD) ratio median {np.nanmedian(r_tail):.3f}")
print("h at n=30: %.3f  c: %.4f" % (0.5 * 30 ** (-1/6), 1 / np.sqrt(1 + (0.5 * 30 ** (-1/6)) ** 2)))
