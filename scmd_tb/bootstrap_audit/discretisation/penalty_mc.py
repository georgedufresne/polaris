"""Bootstrap the 3D fits (K0 uncertainty), Monte Carlo the lattice penalty for a randomly
placed optimum (noise-free and with per-node noise), and a discretisation-aware null."""
import sys, json
import numpy as np
import pandas as pd
sys.path.insert(0, '/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/discretisation')
from common import load, Model, STEP, MU0, E0, AGE0

OUT = '/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/discretisation'
rng = np.random.default_rng(11)
U, V, W = np.meshgrid(np.arange(-8, 10), np.arange(-8, 10), np.arange(-2, 4), indexing='ij')
NODES = np.column_stack([U.ravel(), V.ravel(), W.ravel()]).astype(float)
NMC = 4000
X0R = rng.uniform(0, 1, size=(NMC, 3))
null_split = pd.read_csv('/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02_split/bootstrap_2d_ks_results.csv').ks_statistic.to_numpy()
null_shared = pd.read_csv('/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02/bootstrap_2d_ks_results.csv').ks_statistic.to_numpy()


def qall(A, x0s):
    d = NODES[None, :, :] - x0s[:, None, :]
    return np.einsum('nij,jk,nik->ni', d, A, d)          # (n_x0, n_nodes)


def ks_from_q(kind, K0, p, q):
    if kind == 'cone':
        return K0 + np.sqrt(q)
    if kind == 'rcone':
        return K0 - p + np.sqrt(p ** 2 + q)
    return (K0 ** p + q) ** (1 / p)


d = load()
windows = {'W1 [10.8-12.1, KS<0.04]': ((10.8, 12.1), 0.04), 'W2 [10.8-12.1, KS<0.05]': ((10.8, 12.1), 0.05),
           'W3 [10.4-12.6, KS<0.05]': ((10.4, 12.6), 0.05), 'W4 [11.0-11.8, KS<0.04]': ((11.0, 11.8), 0.04)}
NB = 300
summary = {}
for wname, ((a0, a1), T) in windows.items():
    s = d[(d.age >= a0 - 1e-9) & (d.age <= a1 + 1e-9) & (d.ks < T)].reset_index(drop=True)
    X = s[['u', 'v', 'w']].to_numpy(float)
    y = s.ks.to_numpy()
    for kind in ['hyp', 'rcone']:
        m = Model(kind, 3)
        best = m.fit(X, y, n_starts=24, seed=3)
        K0, x0, A, p = m.unpack(best.x); K0 = abs(K0)
        node0 = m.predict(best.x, np.zeros((1, 3)))[0]
        # lattice: best predicted node over the whole lattice (scored or not)
        q0 = qall(A, x0[None, :])[0]
        lat_best = ks_from_q(kind, K0, p, q0.min())
        # bootstrap
        bs = []
        for b in range(NB):
            idx = rng.integers(0, len(y), len(y))
            r = m.fit(X[idx], y[idx], start=best.x, n_starts=3, seed=100 + b)
            if r is None:
                continue
            k, xx, AA, pp = m.unpack(r.x)
            bs.append((abs(k), AGE0 + xx[0] * STEP[0], MU0 + xx[1] * STEP[1], E0 + xx[2] * STEP[2],
                       m.predict(r.x, np.zeros((1, 3)))[0],
                       ks_from_q(kind, abs(k), pp, qall(AA, xx[None, :])[0].min())))
        bs = np.array(bs)
        # random placement, noise-free and noisy (sigma per node), all nodes assumed scored
        q = qall(A, X0R)
        ksn = ks_from_q(kind, K0, p, q)
        pen0 = ksn.min(axis=1) - K0
        noisy = {}
        for sig in (0.0005, 0.001):
            z = rng.normal(0, sig, size=ksn.shape)
            noisy[sig] = (ksn + z).min(axis=1) - K0
        res = dict(n=len(y), rms=float(np.sqrt(np.mean(best.fun ** 2))), K0=float(K0), eps_or_p=float(p),
                   x0=[float(AGE0 + x0[0] * STEP[0]), float(MU0 + x0[1] * STEP[1]), float(E0 + x0[2] * STEP[2])],
                   pred_best_scored_node=float(node0), lattice_best_pred=float(lat_best),
                   realised_lattice_pen=float(lat_best - K0),
                   K0_boot_median=float(np.median(bs[:, 0])), K0_boot_std=float(np.std(bs[:, 0])),
                   K0_boot_16_84=[float(np.percentile(bs[:, 0], 16)), float(np.percentile(bs[:, 0], 84))],
                   age0_boot_16_84=[float(np.percentile(bs[:, 1], 16)), float(np.percentile(bs[:, 1], 84))],
                   E0_boot_16_84=[float(np.percentile(bs[:, 3], 16)), float(np.percentile(bs[:, 3], 84))],
                   node_pred_boot_16_84=[float(np.percentile(bs[:, 4], 16)), float(np.percentile(bs[:, 4], 84))],
                   lattice_pen_boot_16_84=[float(np.percentile(bs[:, 5] - bs[:, 0], 16)), float(np.percentile(bs[:, 5] - bs[:, 0], 84))],
                   random_pen_mean=float(pen0.mean()), random_pen_median=float(np.median(pen0)),
                   random_pen_p90=float(np.percentile(pen0, 90)), random_pen_max=float(pen0.max()),
                   noisy_gridmin_minus_K0_mean={str(k): float(v.mean()) for k, v in noisy.items()},
                   sqrt_eig=np.sqrt(np.clip(np.linalg.eigvalsh(A), 0, None)).tolist())
        # discretisation-aware nulls (hyperbolic combination of a null draw and the offset term)
        if kind == 'hyp':
            for nm, nul in (('split', null_split), ('shared', null_shared)):
                b0 = nul[rng.integers(0, len(nul), size=q.shape)]
                one = np.sqrt(b0[:, 0] ** 2 + q.min(axis=1) * 0)        # placeholder for the standard null
                # (a) one node: the nearest-in-q node only (no selection among nodes)
                iq = q.argmin(axis=1)
                ka = np.sqrt(b0[np.arange(len(q)), iq] ** 2 + q[np.arange(len(q)), iq])
                # (b) min over all nodes with independent draws (selection included)
                kb = np.sqrt(b0 ** 2 + q).min(axis=1)
                res[f'null_{nm}'] = dict(standard=[float(nul.mean()), float(nul.std()), float(nul.mean() + 3 * nul.std())],
                                         offnode_one_node=[float(ka.mean()), float(ka.std()), float(ka.mean() + 3 * ka.std()), float(np.percentile(ka, 99.865))],
                                         offnode_min_over_nodes=[float(kb.mean()), float(kb.std()), float(kb.mean() + 3 * kb.std()), float(np.percentile(kb, 99.865))])
        summary[f'{wname} {kind}'] = res
        print(f'{wname} {kind}:', json.dumps(res), flush=True)
json.dump(summary, open(f'{OUT}/penalty_mc.json', 'w'), indent=1)
