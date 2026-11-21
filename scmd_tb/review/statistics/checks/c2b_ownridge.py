"""Check 2b: as c2 but each model sample cut around its own traced ridge (as make_scmd_ast does)."""
import sys
exec(open("c2_radial.py").read().split("ridge = ctx.rcut.ridge")[0])
rc = ctx.rcut
for name, model in [("corrected", ast), ("observed", observed_model)]:
    vt, it, rad, ok, vo, io = draw(np.random.RandomState(11), 1_500_000, model)
    tb = ok & rc.in_trace_box(vo - io, vo)
    first = np.flatnonzero(tb)[:rc.n_trace]
    own = rc.trace((vo - io)[first], vo[first])
    sel = ok & rc.select(vo - io, vo, own)
    so, sm = np.ones(len(v), bool), sel & (vo >= 15.02) & (vo < 19.02)
    fr = lambda rr, cut: np.mean(rr < cut)
    n_o = len(R)
    print(f"own ridge, radial={name}: " + "  ".join(
        f"r<{c}: obs {fr(R, c):.4f} model {fr(rad[sm], c):.4f} ({fr(rad[sm], c)/fr(R, c)-1:+.3f}, {(fr(rad[sm], c)-fr(R, c))/np.sqrt(fr(R, c)*(1-fr(R, c))/n_o):+.1f} sigma)"
        for c in (200, 400, 800)))
