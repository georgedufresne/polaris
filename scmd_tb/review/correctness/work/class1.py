import sys; sys.path.insert(0, ".")
from common import *
S = load("dev")
from isochronetoolbox.scmd import ast as AM
s = setup_fresh(S); ctx = S.context(s); ne = ctx.ast.nearest
captured = []
orig = AM.NearestErrors.draw
def spy(self, rs, radius, v, i):
    out = orig(self, rs, radius, v, i)
    captured.append((np.array(radius), np.array(v), np.array(i), out[0]))
    return out
AM.NearestErrors.draw = spy
_, make = ctx.model(10099, 11.4)
cat = make(np.random.RandomState(0), 13.38, 0.015)
AM.NearestErrors.draw = orig
rad, v, i, ok = (np.concatenate([c[k] for c in captured]) for k in range(4))
cls = 2 * (v < ne.split[0]) + (i < ne.split[1])
print("synthetic stars by class:", {int(k): int((cls == k).sum()) for k in np.unique(cls)})
m = cls == 1
if m.any():
    V, I = v[m] + s.zp_v, i[m] + s.zp_i
    print("class-1 synthetic stars: V %.2f-%.2f, I %.2f-%.2f, V-I %.2f-%.2f" % (V.min(), V.max(), I.min(), I.max(), (V-I).min(), (V-I).max()))
    _, near = ne.tree.query(ne._coords(v[m], i[m], rad[m]))
    a = ctx.ast._full
    AV, AI = a["InputF606W"].values[near] + s.zp_v, a["InputF814W"].values[near] + s.zp_i
    print("  their nearest AST: V %.2f-%.2f (median %.2f), I %.2f-%.2f (median %.2f)" % (AV.min(), AV.max(), np.median(AV), AI.min(), AI.max(), np.median(AI)))
    print("  recovered fraction: class-1 %.3f vs all others %.3f" % (ok[m].mean(), ok[~m].mean()))
    rec_nb = ne.rec[ne.nbr[near]].mean(axis=1)
    print("  completeness of the neighbourhoods they get: median %.3f" % np.median(rec_nb))
    # where do class-1 stars sit relative to the window and the sCMD selection?
    Vo = V
    print("  class-1 stars inside the V window 15.02-19.02:", int(((Vo > 15.02) & (Vo < 19.02)).sum()))
