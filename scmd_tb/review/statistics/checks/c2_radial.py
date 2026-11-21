"""Check 2: completeness weights and the corrected radial distribution, end to end, on the newmodel setup."""
import sys, json, time, numpy as np, pandas as pd
sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/review/statistics/dev")
from scipy.interpolate import interp1d
from isochronetoolbox.scmd import Setup, context
from isochronetoolbox.scmd.ast import ArtificialStars, centre
from isochronetoolbox.scmd.model import select_age, apparent
F32 = np.float32
S = Setup(**{k: v for k, v in json.load(open("/dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/newmodel/setup.json")).items()
             if k in Setup.__dataclass_fields__})
t = time.time(); ctx = context(S); print(f"context {time.time()-t:.0f} s", flush=True)
zpv, zpi = S.zp_v, S.zp_i
ast = ctx.ast
# --- observed fit stars (the ridge selection the context made) with positions
raw = pd.read_csv(S.obs_raw, sep=r"\s+")
keep = ctx.rcut.select(raw.VIvega.values, raw.Vvega.values, ctx.rcut.ridge)
d = raw[keep]
v, i = d.Vvega.values.astype(float), d.Ivega.values.astype(float)
x, y = d.x.values, d.y.values
fxc, fyc = centre(x.astype(F32), y.astype(F32))
R = np.hypot(x.astype(F32) - F32(fxc), y.astype(F32) - F32(fyc)).astype(float)
A = ast._full; axc, ayc = centre(A.InputX.values, A.InputY.values)
print(f"fit stars {len(v)}; fit-star centre ({fxc:.1f},{fyc:.1f}) AST centre ({axc:.1f},{ayc:.1f}) offset {np.hypot(fxc-axc, fyc-ayc):.1f} px")
ne = ast.nearest
c100 = ne.completeness(v - zpv, i - zpi, R)
_, near = ne.tree.query(ne._coords(v - zpv, i - zpi, R))
c30 = ne.rec[ne.nbr[near]].mean(1)
print("c100 at fit stars pct 0.1,1,5,50:", np.round(np.percentile(c100, [0.1, 1, 5, 50]), 3), " capped (c<0.1):", (c100 < 0.1).sum())
for lo, hi in [(0, 100), (100, 200), (200, 400), (400, 800), (800, 1500), (1500, 5000)]:
    s = (R >= lo) & (R < hi)
    print(f"  r {lo:4d}-{hi:4d}: n {s.sum():6d} mean c100 {c100[s].mean():.3f} mean c30(pick) {c30[s].mean():.3f}  E[1/c100] {np.mean(1/np.maximum(c100[s],.1)):.3f}")
# --- model stars for the best cell, radius kept
run, age, dm, ebv = 10099, 11.4, 13.38, 0.015
age_yr, rows = select_age(ctx.isochrones(run), age, S.exact_age)
ph = apparent(ctx.absolute(rows), dm, ebv, {"F606W": 1.0, "F814W": 0.0}, "F606W")
mass = rows["initial_mass"].values; o = np.argsort(mass)
m, Vs, Is = mass[o], ph["F606W"][o], ph["F814W"][o]
lo, hi = ctx.v_window
w = np.where((Vs >= lo - S.pad) & (Vs <= hi + S.pad))[0]; a_, b_ = w[0], w[-1]
m_lo, m_hi = m[a_:b_ + 1].min(), m[a_:b_ + 1].max()
fV = interp1d(m, Vs, kind="linear", fill_value="extrapolate"); fI = interp1d(m, Is, kind="linear", fill_value="extrapolate")
e = S.imf_slope + 1.0
def draw(rs, n, model):
    u = rs.random_sample(n); ms = (u * (m_hi ** e - m_lo ** e) + m_lo ** e) ** (1 / e)
    vt, it = fV(ms), fI(ms)
    bb = rs.random_sample(n) <= S.binary_fraction
    m2 = (0.5 + 0.5 * rs.random_sample(int(bb.sum()))) * ms[bb]
    vt[bb] = -2.5 * np.log10(10 ** (-0.4 * vt[bb]) + 10 ** (-0.4 * fV(m2)))
    it[bb] = -2.5 * np.log10(10 ** (-0.4 * it[bb]) + 10 ** (-0.4 * fI(m2)))
    rad = model.radii(rs, n, vt)
    ok, ev, ei = model.nearest.draw(rs, rad, vt - zpv, it - zpi)
    return vt, it, rad, ok, vt - ev, it - ei
observed_model = ArtificialStars(ast._full, *ast._args, k=S.ast_k, smooth=True, metric=S.ast_metric, radial="observed")
ridge = ctx.rcut.ridge
for name, model in [("corrected", ast), ("observed", observed_model)]:
    vt, it, rad, ok, vo, io = draw(np.random.RandomState(11), 1_500_000, model)
    sel = ok & ctx.rcut.select(vo - io, vo, ridge)
    print(f"--- radial={name}: model selected {sel.sum()}")
    # P(ridge-select | recovered, true V in window) vs radius
    base = ok & (vt > lo + 0.1) & (vt < hi - 0.1)
    line = []
    for r0, r1 in [(0, 100), (100, 200), (200, 400), (400, 800), (800, 1500), (1500, 5000)]:
        s = base & (rad >= r0) & (rad < r1)
        line.append(f"{r0}-{r1}: {np.mean(sel[s]):.3f}")
    print("  s(r) = P(pass ridge | recovered, true V inside window): " + "  ".join(line))
    for b0, b1 in [(15.02, 16.0), (16.0, 17.0), (17.0, 18.0), (18.0, 19.02), (15.02, 19.02)]:
        so = (v >= b0) & (v < b1); sm = sel & (vo >= b0) & (vo < b1)
        fr = lambda rr, cut: np.mean(rr < cut)
        print(f"  V {b0:5.2f}-{b1:5.2f}: frac r<200 obs {fr(R[so],200):.4f} model {fr(rad[sm],200):.4f} | r<400 obs {fr(R[so],400):.4f} "
              f"model {fr(rad[sm],400):.4f} | r<800 obs {fr(R[so],800):.4f} model {fr(rad[sm],800):.4f} | median r obs {np.median(R[so]):.0f} model {np.median(rad[sm]):.0f}")
