"""Round 2: is the table's split (seed 0) unusual, and more splits / AST resamples.
  P1  split seed 0 with fresh refs (seeds 1000, 1005): half-B0 reals vs A0 refs (split) and F refs (rmirror)
  P2  full-list reals vs A0/F refs (seeds 1000, 1005) and A_s/F_s for s=5..8 (mirror, shared)
  P3  half-B_s reals for s=5,6 vs A_s (split), F_s (rmirror); AST-bootstrap lists b=5,6 vs F_b (astboot)
  P4  same for s=7,8
"""
import copy, os, sys
from common import *
role, K = sys.argv[1], int(sys.argv[2])
out_csv = os.path.join(OUT, f'results2_{role}.csv')
T0 = time.time()
rows = []
def log(m):
    print(f'[{role} {time.time() - T0:7.1f}s] {m}', flush=True)
def save(design, tag, seed, ks, n_sim, n_ref):
    rows.append(dict(role=role, design=design, tag=str(tag), seed=int(seed), ks=ks, n_sim=n_sim, n_ref=n_ref))
    pd.DataFrame(rows).to_csv(out_csv, index=False)
def resampled_copy(ast, ast_csv, b):
    a = pd.read_csv(ast_csv, dtype=np.float32)
    idx = np.random.RandomState(10_000 + b).randint(0, len(a), len(a))
    out = copy.copy(ast)
    out.nearest = NearestErrors(a.iloc[idx].reset_index(drop=True))
    return out
S = setup(); ctx = context(S); full = ctx.ast; _, make = ctx.model(RUN, AGE)
log('context ready')
if role in ('P1', 'P2'):
    ha0, hb0 = full.halves(0)
    refs = [('0r1000', 'A', Ref(sample(ctx, make, ha0, 1000), 'A0_1000')), ('0r1000', 'F', Ref(sample(ctx, make, full, 1000), 'F_1000')),
            ('0r1005', 'A', Ref(sample(ctx, make, ha0, 1005), 'A0_1005')), ('0r1005', 'F', Ref(sample(ctx, make, full, 1005), 'F_1005'))]
    if role == 'P1':
        for i in range(K):
            sim = sample(ctx, make, hb0, i)
            ks = ks_multi(sim, [r for _, _, r in refs])
            for (tag, d, r), k in zip(refs, ks):
                save('split' if d == 'A' else 'rmirror', tag, i, k, len(sim), len(r.pts))
            log(f'i={i}: ' + ' '.join(f'{t}{d} {k:.4f}' for (t, d, _), k in zip(refs, ks)))
    else:
        for s in (5, 6, 7, 8):
            ha, _ = full.halves(s)
            refs += [(str(s), 'A', Ref(sample(ctx, make, ha, 1000 + s), f'A{s}')), (str(s), 'F', Ref(sample(ctx, make, full, 1000 + s), f'F{s}'))]
        log('refs built')
        for i in range(K):
            sim = sample(ctx, make, full, i)
            ks = ks_multi(sim, [r for _, _, r in refs])
            for (tag, d, r), k in zip(refs, ks):
                save('mirror' if d == 'A' else 'shared', tag, i, k, len(sim), len(r.pts))
            log(f'i={i}: ' + ' '.join(f'{t}{d} {k:.4f}' for (t, d, _), k in zip(refs, ks)))
else:
    ss = (5, 6) if role == 'P3' else (7, 8)
    for s in ss:
        ha, hb = full.halves(s)
        rA, rF = Ref(sample(ctx, make, ha, 1000 + s), f'A{s}'), Ref(sample(ctx, make, full, 1000 + s), f'F{s}')
        for i in range(K):
            sim = sample(ctx, make, hb, i)
            k1, k2 = ks_multi(sim, [rA, rF])
            save('split', s, i, k1, len(sim), len(rA.pts)); save('rmirror', s, i, k2, len(sim), len(rF.pts))
            log(f's={s} i={i}: split {k1:.5f} rmirror {k2:.5f}')
        lst = resampled_copy(full, S.ast, s)
        for i in range(K):
            sim = sample(ctx, make, lst, i)
            k, = ks_multi(sim, [rF])
            save('astboot', s, i, k, len(sim), len(rF.pts))
            log(f'b={s} i={i}: astboot {k:.5f}')
log('done')
