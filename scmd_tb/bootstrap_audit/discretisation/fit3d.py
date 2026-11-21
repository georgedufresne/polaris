"""3D fits of KS over (age, mu, E(V-I)) near the minimum for run 10099, using the
all-cells table at 11.4 Gyr plus the stage-3 rows (best four cells by JS) at other ages.
Model families: hyp  KS^2 = K0^2 + q;  cone KS = K0 + sqrt(q);  rcone KS = K0 - eps + sqrt(eps^2 + q);
par KS = K0 + q;  pfree KS^p = K0^p + q.  q = (x - x0)^T A (x - x0), x in grid steps."""
import sys, json
import numpy as np
sys.path.insert(0, '/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/discretisation')
from common import load, Model, STEP, MU0, E0, AGE0

rng = np.random.default_rng(1)
# nodes around the origin for lattice minimisation (grid steps)
U, V, W = np.meshgrid(np.arange(-12, 13), np.arange(-12, 13), np.arange(-3, 4), indexing='ij')
NODES = np.column_stack([U.ravel(), V.ravel(), W.ravel()]).astype(float)
X0R = rng.uniform(0, 1, size=(20000, 3))          # random optimum positions within one cell


def qmin_lattice(A, x0s):
    out = np.empty(len(x0s))
    for i in range(0, len(x0s), 500):
        d = NODES[None, :, :] - x0s[i:i + 500, None, :]
        q = np.einsum('nij,jk,nik->ni', d, A, d)
        out[i:i + 500] = q.min(axis=1)
    return out


def penalty(m, par, q):
    K0, x0, A, p = m.unpack(par)
    K0 = abs(K0)
    if m.kind == 'cone':
        return np.sqrt(q)
    if m.kind == 'rcone':
        return -p + np.sqrt(p ** 2 + q)
    return (K0 ** p + q) ** (1 / p) - K0


def describe(m, r, X, y, label):
    K0, x0, A, p = m.unpack(r.x)
    K0 = abs(K0)
    rms = np.sqrt(np.mean(r.fun ** 2))
    ev, evec = np.linalg.eigh(A)
    # realised: model KS at the best node of the lattice (noise-free) vs K0
    pred_nodes = m.predict(r.x, NODES)
    ib = np.argmin(pred_nodes)
    best_node = NODES[ib]
    q_r = qmin_lattice(A, x0[None, :])[0]
    pen_r = penalty(m, r.x, np.array([q_r]))[0]
    pen = penalty(m, r.x, qmin_lattice(A, X0R))
    pred0 = m.predict(r.x, np.zeros((1, 3)))[0]
    phys = AGE0 + x0[0] * STEP[0], MU0 + x0[1] * STEP[1], E0 + x0[2] * STEP[2]
    pstr = '' if p is None else f' p/eps={p:.4g}'
    print(f'{label} {m.kind:5s} n={len(y):3d} rms={rms:.5f} K0={K0:.5f}{pstr} x0=(age {phys[0]:.3f}, mu {phys[1]:.4f}, E {phys[2]:.5f}) '
          f'pred@node={pred0:.5f} best_node={best_node.astype(int).tolist()} realised_pen={pen_r:.5f} '
          f'random_pen mean={pen.mean():.5f} med={np.median(pen):.5f} p90={np.percentile(pen, 90):.5f} '
          f'sqrt(eig A)/step={np.sqrt(np.clip(ev, 0, None)).round(5).tolist()}')
    return dict(label=label, kind=m.kind, n=int(len(y)), rms=float(rms), K0=float(K0), p=None if p is None else float(p),
                age0=float(phys[0]), mu0=float(phys[1]), E0=float(phys[2]), pred_node=float(pred0),
                best_node=best_node.astype(int).tolist(), realised_pen=float(pen_r),
                pen_mean=float(pen.mean()), pen_med=float(np.median(pen)), pen_p90=float(np.percentile(pen, 90)),
                sqrt_eig=np.sqrt(np.clip(ev, 0, None)).tolist(), evec=evec.tolist(), par=r.x.tolist())


if __name__ == '__main__':
    d = load()
    out = []
    windows = [((11.0, 11.8), 0.04), ((10.8, 12.1), 0.04), ((10.8, 12.1), 0.05), ((10.8, 12.1), 0.06),
               ((10.4, 12.6), 0.05), ((11.0, 11.8), 0.06)]
    kinds = sys.argv[1].split(',') if len(sys.argv) > 1 else ['hyp', 'cone', 'rcone', 'par', 'pfree']
    for (a0, a1), T in windows:
        s = d[(d.age >= a0 - 1e-9) & (d.age <= a1 + 1e-9) & (d.ks < T)]
        X = s[['u', 'v', 'w']].to_numpy(float)
        y = s.ks.to_numpy()
        label = f'[{a0}-{a1} Gyr, KS<{T}]'
        for kind in kinds:
            m = Model(kind, k=3)
            r = m.fit(X, y, n_starts=24, seed=3)
            out.append(describe(m, r, X, y, label))
        sys.stdout.flush()
    json.dump(out, open('/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/discretisation/fit3d.json', 'w'), indent=1)
