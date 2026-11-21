import sys, time; sys.path.insert(0, ".")
from common import *
S = load("dev")
from isochronetoolbox.scmd import ast as AM
s = setup_fresh(S); ctx = S.context(s); A = ctx.ast
# capture the arguments make_scmd_ast passes to nearest.draw for one real sCMD of the bootstrap cell
captured = []
orig = AM.NearestErrors.draw
def spy(self, rs, radius, v, i):
    captured.append((np.array(radius), np.array(v), np.array(i)))
    return orig(self, rs, radius, v, i)
AM.NearestErrors.draw = spy
_, make = ctx.model(10099, 11.4)
cat = make(np.random.RandomState(0), 13.38, 0.015)
AM.NearestErrors.draw = orig
rad = np.concatenate([c[0] for c in captured]); v = np.concatenate([c[1] for c in captured]); i = np.concatenate([c[2] for c in captured])
print("synthetic stars drawn for one sCMD:", len(rad))
def usage(ne, tag):
    _, near = ne.tree.query(ne._coords(v, i, rad))
    n = ne.rec[ne.nbr[near]].sum(axis=1)
    rs = np.random.RandomState(1)
    ok, ev, ei = ne.draw(rs, rad, v, i)
    nok = n[ok]
    sv = ne._L[near, 0][ok] / (ne.bw * np.maximum(nok, 1) ** (-1/6))
    print(f"{tag}: recovered synthetic stars {ok.sum()}; in neighbourhoods with n=1 (zero-width kernel): {(nok==1).sum()} "
          f"({100*(nok==1).mean():.2f}%), n<=4: {100*(nok<=4).mean():.2f}%, MAD floor: {100*np.isclose(sv,1e-4).mean():.2f}%")
    return ok, ev, ei
ok_s, ev_s, ei_s = usage(A.nearest, "full AST, smooth")
R = A.resample(1_000_007)
usage(R.nearest, "resample, smooth")
# copy vs smooth on the same synthetic stars and the same pick stream
cp = AM.NearestErrors.__new__(AM.NearestErrors); cp.__dict__.update(A.nearest.__dict__); cp.smooth = False
ok_c, ev_c, ei_c = cp.draw(np.random.RandomState(1), rad, v, i)
assert np.array_equal(ok_c, ok_s)
dc, ds = ev_c[ok_c], ev_s[ok_s]
mad = lambda x: 1.4826 * np.median(np.abs(x - np.median(x)))
print("dV copy vs smooth: MAD %.4f / %.4f; SD %.4f / %.4f; mean %.5f / %.5f" % (mad(dc), mad(ds), dc.std(), ds.std(), dc.mean(), ds.mean()))
for q in (0.1, 0.2, 0.3):
    print(f"  P(dV > {q}) copy {np.mean(dc > q):.5f} smooth {np.mean(ds > q):.5f}   P(dV < -{q}) copy {np.mean(dc < -q):.5f} smooth {np.mean(ds < -q):.5f}")
print("  0.1/99.9 percentiles copy", np.percentile(dc, [0.1, 99.9]).round(4), "smooth", np.percentile(ds, [0.1, 99.9]).round(4))
col_c, col_s = (ev_c - ei_c)[ok_c], (ev_s - ei_s)[ok_s]
print("colour error MAD copy %.4f smooth %.4f; SD copy %.4f smooth %.4f" % (mad(col_c), mad(col_s), col_c.std(), col_s.std()))
# repeated values in smooth draws (the docstring: 'no two sCMDs share error values')
u, cnt = np.unique(ds, return_counts=True)
print("smooth dV values that repeat within one batch:", int((cnt > 1).sum()), "values,", int(cnt[cnt > 1].sum()), "stars")
