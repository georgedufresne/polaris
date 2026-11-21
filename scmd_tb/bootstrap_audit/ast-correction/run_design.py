"""ast-correction audit: direct simulations of the bootstrap null with different AST-list designs.

Roles (one process each):
  V  verify the exact bootstrap code path (bootstrap.realisations, split and shared) against the
     10,000-realisation tables, capture the two table references, then score extra realisations:
       full-list realisations vs table-shared ref (reproduces shared table) and vs table-split ref (mirror)
       half-B(seed 0) realisations vs table-split ref (reproduces split table)
  F  full-list realisations vs fresh references, splits s=1..4:
       ref_A_s (half A_s)      -> 'mirror'   (AST term ~ (e_A - e_B)/2, variance = full-list AST noise)
       ref_F_s (full list)     -> 'shared'
       ref_S_s (full, smoothed)-> 'smoothed' (task's suggested smoothed error model), s=1,2
       ref_SA_s (half A_s, smoothed) -> 'smoothedA' (same smoothing, independent list), s=1,2
  B  half-B_s realisations vs ref_A_s ('split') and vs ref_F_s ('rmirror': full ref vs half-B, also ~ (e_A-e_B)/2)
  C  realisations from a bootstrap-resampled AST list L*_b vs ref_F_b ('astboot')
"""
import copy, os, sys
from common import *
import isochronetoolbox.scmd.bootstrap as BOOT

role, K = sys.argv[1], int(sys.argv[2])
SPLITS = [1, 2, 3, 4]
REFSEED = 1000          # fresh reference seeds 1000 + s (realisations use 0..K-1, as the tables)
H = 0.5                 # smoothing bandwidth / local robust scale
out_csv = os.path.join(OUT, f'results_{role}.csv')
T0 = time.time()


def log(m):
    print(f'[{role} {time.time() - T0:7.1f}s] {m}', flush=True)


rows = []


def save(design, tag, seed, ks, n_sim, n_ref):
    rows.append(dict(role=role, design=design, tag=tag, seed=int(seed), ks=ks, n_sim=n_sim, n_ref=n_ref))
    pd.DataFrame(rows).to_csv(out_csv, index=False)


def smoothed_copy(ast, h=H, jitter_seed=0):
    """A copy of ``ast`` whose NearestErrors.draw adds variance-preserving jitter: pick as usual,
    then dV,dI -> med + c (pick - med) + c h L z, c = 1/sqrt(1+h^2), med and L (robust, MAD-based,
    with the dV-dI correlation) from the recovered errors of the nearest AST's 30 neighbours.
    The jitter uses its own RandomState so the star sampling matches the unsmoothed reference."""
    ne0 = ast.nearest
    ne = copy.copy(ne0)
    nbr = ne0.nbr
    okn = ne0.rec[nbr]
    EV = np.where(okn, ne0.ev[nbr], np.nan)
    EI = np.where(okn, ne0.ei[nbr], np.nan)
    mv, mi = np.nanmedian(EV, 1), np.nanmedian(EI, 1)
    sv = 1.4826 * np.nanmedian(np.abs(EV - mv[:, None]), 1)
    si = 1.4826 * np.nanmedian(np.abs(EI - mi[:, None]), 1)
    sc = 1.4826 * np.nanmedian(np.abs((EV - EI) - (mv - mi)[:, None]), 1)
    rho = np.clip((sv ** 2 + si ** 2 - sc ** 2) / np.maximum(2 * sv * si, 1e-12), -0.99, 0.99)
    bad = (okn.sum(1) < 5) | ~np.isfinite(sv) | ~np.isfinite(si)
    sv[bad] = si[bad] = 0.0
    mv[bad] = mi[bad] = 0.0
    rho[bad] = 0.0
    jr = np.random.RandomState(jitter_seed)
    c = 1.0 / np.sqrt(1 + h * h)

    def draw(self, rs, radius, v_inst, i_inst):
        _, near = self.tree.query(np.column_stack([v_inst, i_inst, radius]) / self.scale)
        pick = self.nbr[near, np.minimum((rs.random_sample(len(near)) * self.k).astype(int), self.k - 1)]
        z1, z2 = jr.standard_normal(len(near)), jr.standard_normal(len(near))
        b = bad[near]
        dv = np.where(b, self.ev[pick], mv[near] + c * (self.ev[pick] - mv[near]) + c * h * sv[near] * z1)
        di = np.where(b, self.ei[pick], mi[near] + c * (self.ei[pick] - mi[near])
                      + c * h * si[near] * (rho[near] * z1 + np.sqrt(1 - rho[near] ** 2) * z2))
        return self.rec[pick], dv, di
    ne.draw = types.MethodType(draw, ne)
    out = copy.copy(ast)
    out.nearest = ne
    return out


def resampled_copy(ast, ast_csv, b):
    """``ast`` with its NearestErrors rebuilt from a bootstrap resample (with replacement, same size)
    of the artificial-star list: a synthetic AST list drawn from the full list as 'truth'."""
    a = pd.read_csv(ast_csv, dtype=np.float32)
    idx = np.random.RandomState(10_000 + b).randint(0, len(a), len(a))
    out = copy.copy(ast)
    out.nearest = NearestErrors(a.iloc[idx].reset_index(drop=True))
    return out, len(np.unique(idx))


S = setup()
ctx = context(S)
full = ctx.ast
_, make = ctx.model(RUN, AGE)
log(f'context ready, n_obs {len(ctx.obs3)}')

if role == 'V':
    sh = pd.read_csv(G + '/bootstrap_2d_ks_results.csv').set_index('bootstrap_number')['ks_statistic']
    sp = pd.read_csv(G + '_split/bootstrap_2d_ks_results.csv').set_index('bootstrap_number')['ks_statistic']
    cap = {}
    orig = BOOT.reference

    def capture(c, run, age, dm, ebv):
        r = orig(c, run, age, dm, ebv)
        cap['ref'], cap['n_ast'] = r, c.ast.n_ast
        return r
    BOOT.reference = capture
    try:
        for split, table, name in [(True, sp, 'split'), (False, sh, 'shared')]:
            before = ctx.ast
            rr = BOOT.realisations(ctx, RUN, AGE, DM, EVI, ids=[0, 1, 2], split_ast=split)
            log(f'bootstrap.realisations split_ast={split}: reference built with n_ast={cap["n_ast"]} '
                f'(full {full.n_ast}); ctx.ast restored: {ctx.ast is before}')
            for r in rr:
                i = r['bootstrap_number']
                log(f'   {name} row {i}: code {r["ks_statistic"]:.10f}  table {table[i]:.10f}  diff {r["ks_statistic"] - table[i]:.2e}')
                save(f'verify_{name}', 'table', i, r['ks_statistic'], r['n_simulated'], r['n_observed'])
            cap[name] = Ref(cap['ref'], name)
    finally:
        BOOT.reference = orig
    refT_split, refT_shared = cap['split'], cap['shared']
    ha0, hb0 = full.halves(0)
    for i in range(K):
        sim = sample(ctx, make, full, i)
        k_sh, k_mir = ks_multi(sim, [refT_shared, refT_split])
        save('shared', 'T', i, k_sh, len(sim), len(refT_shared.pts))
        save('mirror', 'T', i, k_mir, len(sim), len(refT_split.pts))
        simb = sample(ctx, make, hb0, i)
        k_sp, = ks_multi(simb, [refT_split])
        save('split', 'T', i, k_sp, len(simb), len(refT_split.pts))
        log(f'i={i}: shared {k_sh:.5f} (table {sh[i]:.5f})  split {k_sp:.5f} (table {sp[i]:.5f})  mirror {k_mir:.5f}')

elif role in ('F', 'B'):
    refs = {}
    for s in SPLITS:
        ha, hb = full.halves(s)
        refs[s] = dict(A=Ref(sample(ctx, make, ha, REFSEED + s), f'A{s}'),
                       F=Ref(sample(ctx, make, full, REFSEED + s), f'F{s}'), hb=hb)
        if role == 'F' and s in (1, 2):
            refs[s]['S'] = Ref(sample(ctx, make, smoothed_copy(full, H, 5000 + s), REFSEED + s), f'S{s}')
            refs[s]['SA'] = Ref(sample(ctx, make, smoothed_copy(ha, H, 5000 + s), REFSEED + s), f'SA{s}')
    log('references built')
    if role == 'F':
        order = [(s, d) for s in SPLITS for d in ('A', 'F', 'S', 'SA') if d in refs[s]]
        name = dict(A='mirror', F='shared', S='smoothed', SA='smoothedA')
        for i in range(K):
            sim = sample(ctx, make, full, i)
            ks = ks_multi(sim, [refs[s][d] for s, d in order])
            for (s, d), k in zip(order, ks):
                save(name[d], s, i, k, len(sim), len(refs[s][d].pts))
            log(f'i={i}: ' + '  '.join(f'{name[d][:3]}{s} {k:.4f}' for (s, d), k in zip(order, ks)))
    else:
        for s in SPLITS:
            for i in range(K):
                sim = sample(ctx, make, refs[s]['hb'], i)
                k_sp, k_rm = ks_multi(sim, [refs[s]['A'], refs[s]['F']])
                save('split', s, i, k_sp, len(sim), len(refs[s]['A'].pts))
                save('rmirror', s, i, k_rm, len(sim), len(refs[s]['F'].pts))
                log(f's={s} i={i}: split {k_sp:.5f}  rmirror {k_rm:.5f}')

elif role == 'C':
    for b in SPLITS:
        refF = Ref(sample(ctx, make, full, REFSEED + b), f'F{b}')
        lst, nu = resampled_copy(full, S.ast, b)
        log(f'resample {b}: {nu} unique of {full.n_ast}')
        for i in range(K):
            sim = sample(ctx, make, lst, i)
            k, = ks_multi(sim, [refF])
            save('astboot', b, i, k, len(sim), len(refF.pts))
            log(f'b={b} i={i}: astboot {k:.5f}')
log('done')
