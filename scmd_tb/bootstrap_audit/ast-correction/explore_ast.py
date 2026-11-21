from common import *
t0 = time.time()
S = setup()
ctx = context(S)
print('context %.1f s' % (time.time() - t0), 'obs3', ctx.obs3.shape, 'n_target', ctx.n_target)
full = ctx.ast
ne = full.nearest
a = pd.read_csv(S.ast, dtype=np.float32)
print('AST rows', len(a), 'recovered', ne.rec.sum(), '(%.3f)' % ne.rec.mean())
vin, iin = a.InputF606W.values.astype(float), a.InputF814W.values.astype(float)
V, I = vin + S.zp_v, iin + S.zp_i
print('input V range %.2f..%.2f, I range %.2f..%.2f' % (V.min(), V.max(), I.min(), I.max()))
print('corr(Vin, Iin) %.5f; Vin-Iin spread (std) %.4f' % (np.corrcoef(V, I)[0, 1], np.std(V - I)))
w = (V >= 15.02 - 1) & (V <= 19.02 + 1)
print('ASTs with V_in within window +- 1 mag:', w.sum(), 'recovered', (w & ne.rec).sum())
w0 = (V >= 15.02) & (V <= 19.02)
print('ASTs with V_in within window:', w0.sum(), 'recovered', (w0 & ne.rec).sum())
x, y = a.InputX.values, a.InputY.values
xc, yc = centre(x, y)
r = np.hypot(x - xc, y - yc)
print('AST centre (%.2f, %.2f); radius range %.0f..%.0f px' % (xc, yc, r.min(), r.max()))
print('observed radial distribution (fit stars): pct 5/50/95 =', np.percentile(full.dist, [5, 50, 95]).round(0))
print('AST radius pct 5/50/95 (window)=', np.percentile(r[w0], [5, 50, 95]).round(0))
# errors: correlation of dV, dI among recovered in window
ev, ei = ne.ev, ne.ei
s = w0 & ne.rec
print('corr(dV,dI) recovered in window %.3f' % np.corrcoef(ev[s], ei[s])[0, 1])
dc = ev - ei
print('robust sigma dV %.4f dI %.4f dColour %.4f (MAD*1.4826)' % tuple(1.4826 * np.median(np.abs(z[s] - np.median(z[s]))) for z in (ev, ei, dc)))
# neighbourhood size in scaled units
d, _ = ne.tree.query(ne.tree.data[s][:5000], k=30)
print('30-NN radius (scaled units) median %.3f -> dV ~ %.3f mag, dr ~ %.0f px' % (np.median(d[:, -1]), np.median(d[:, -1]) * 0.1, np.median(d[:, -1]) * 200))
# halves
ha, hb = full.halves()
print('halves n_ast', ha.n_ast, hb.n_ast, 'sum', ha.n_ast + hb.n_ast, 'of', full.n_ast)
n = len(a)
half = np.random.RandomState(0).random_sample(n) < 0.5
print('disjoint by construction: overlap', np.sum(half & ~half), ' union', np.sum(half | ~half))
for nm, h in [('A', ha), ('B', hb)]:
    m = half if nm == 'A' else ~half
    xh, yh = x[m], y[m]
    c = centre(xh, yh)
    print('half %s: centre (%.2f, %.2f) shift %.2f px; nearest.ev identical to full subset: %s; nbr table shape %s, max index %d' % (
        nm, c[0], c[1], np.hypot(c[0] - xc, c[1] - yc), np.array_equal(h.nearest.ev, ne.ev[m]), h.nearest.nbr.shape, h.nearest.nbr.max()))
    dh, _ = h.nearest.tree.query(h.nearest.tree.data[:5000], k=30)
    print('   30-NN radius (scaled) median %.3f' % np.median(dh[:, -1]))
print('dist identical across halves:', np.array_equal(ha.dist, full.dist), np.array_equal(hb.dist, full.dist))
print('elapsed %.1f s' % (time.time() - t0))
