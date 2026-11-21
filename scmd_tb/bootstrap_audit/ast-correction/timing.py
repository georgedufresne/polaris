from common import *
S = setup(); ctx = context(S)
_, make = ctx.model(RUN, AGE)
t = time.time(); p0 = sample(ctx, make, ctx.ast, 0); t1 = time.time() - t
t = time.time(); r0 = Ref(p0, 'r0'); t2 = time.time() - t
t = time.time(); p1 = sample(ctx, make, ctx.ast, 1); t3 = time.time() - t
t = time.time(); k = ks_multi(p1, [r0]); t4 = time.time() - t
t = time.time(); k2 = ks_2d(r0.pts, p1, r0.cdf)[0]; t5 = time.time() - t
t = time.time(); k3 = ks_multi(p1, [r0, r0, r0, r0]); t6 = time.time() - t
print('sCMD %.2f s, rank2d ref %.2f s, sCMD %.2f s, ks_multi(1) %.2f s, ks_2d %.2f s, ks_multi(4) %.2f s' % (t1, t2, t3, t4, t5, t6))
print('n', len(p0), len(p1), 'KS', k, k2, k3)
t = time.time(); ha, hb = ctx.ast.halves(); print('halves %.2f s' % (time.time() - t))
