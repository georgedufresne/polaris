"""Shared loading and model code for the discretisation audit (existing tables only)."""
import numpy as np
import pandas as pd
from scipy.optimize import least_squares

G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
RUN, AGE0, MU0, E0 = 10099, 11.4, 13.38, 0.015
STEP = np.array([0.1, 0.01, 0.005])          # age (Gyr), mu, E(V-I)


def load(run=RUN):
    """All cells at 11.4 Gyr (all_cells table) + stage-3 rows of `run` at the other ages.
    Coordinates in grid steps relative to the best node: u (age), v (mu), w (E)."""
    c = pd.read_csv(f'{G}/combined_2d_ks_results.csv')
    r = c[c.simulation_id == run][['age_gyr', 'distance_modulus', 'extinction_ebv', 'ks_statistic', 'kde_rank']]
    r = r.rename(columns={'age_gyr': 'age', 'distance_modulus': 'mu', 'extinction_ebv': 'E', 'ks_statistic': 'ks'})
    r['src'] = 'stage3'
    if run == RUN:
        a = pd.read_csv(f'{G}/all_cells_{RUN}_{AGE0}.csv').rename(columns={'evi': 'E'})
        a['age'] = AGE0
        a['src'] = 'allcells'
        a['kde_rank'] = a['js_rank']
        r = r[~np.isclose(r.age, AGE0)]
        r = pd.concat([a[['age', 'mu', 'E', 'ks', 'kde_rank', 'src']], r], ignore_index=True)
    r['u'] = np.round((r.age - AGE0) / STEP[0]).astype(int)
    r['v'] = np.round((r.mu - MU0) / STEP[1]).astype(int)
    r['w'] = np.round((r.E - E0) / STEP[2]).astype(int)
    return r.reset_index(drop=True)


def chol(par, k):
    """Lower-triangular L (k x k) from k(k+1)/2 numbers; positive diagonal via exp."""
    L = np.zeros((k, k))
    idx = np.tril_indices(k)
    L[idx] = par
    L[np.diag_indices(k)] = np.exp(np.clip(np.diag(L), -50, 50))
    return L


class Model:
    """KS^p = K0^p + d^T A d (p = 2: hyperbolic, p = 1: paraboloid), or the cone
    KS = K0 + sqrt(d^T A d).  d = x - x0 in grid steps, A = L L^T."""

    def __init__(self, kind='hyp', k=3, p=None):
        self.kind, self.k, self.p = kind, k, p

    def npar(self):
        return 1 + self.k + self.k * (self.k + 1) // 2 + (1 if self.kind in ('pfree', 'rcone') else 0)

    def unpack(self, par):
        k = self.k
        K0 = par[0]
        x0 = par[1:1 + k]
        L = chol(par[1 + k:1 + k + k * (k + 1) // 2], k)
        p = par[-1] if self.kind in ('pfree', 'rcone') else {'hyp': 2.0, 'par': 1.0, 'cone': None}[self.kind]
        return K0, x0, L @ L.T, p

    def q(self, par, X):
        K0, x0, A, p = self.unpack(par)
        d = X - x0
        return np.einsum('ij,jk,ik->i', d, A, d)

    def predict(self, par, X):
        K0, x0, A, p = self.unpack(par)
        q = self.q(par, X)
        if self.kind == 'cone':
            return np.abs(K0) + np.sqrt(q)
        if self.kind == 'rcone':            # p holds the rounding scale eps >= 0
            return np.abs(K0) - p + np.sqrt(p ** 2 + q)
        return (np.abs(K0) ** p + q) ** (1.0 / p)

    def fit(self, X, y, start=None, n_starts=12, seed=0):
        rng = np.random.default_rng(seed)
        best = None
        k = self.k
        for s in range(n_starts):
            if start is not None and s == 0:
                p0 = np.asarray(start, float)
            else:
                K0 = y.min() * rng.uniform(0.5, 1.0)
                x0 = X[np.argmin(y)] + rng.normal(0, 0.3, k)
                Ld = np.zeros(k * (k + 1) // 2)
                # crude diagonal guess, log-scale, in KS units^(p/2) per step
                scale = {'hyp': np.log(0.01), 'par': np.log(np.sqrt(0.003)), 'cone': np.log(0.01), 'pfree': np.log(0.01), 'rcone': np.log(0.01)}[self.kind]
                Ld[[i * (i + 1) // 2 + i for i in range(k)]] = scale + rng.normal(0, 1.0, k)
                extra = [[2.0]] if self.kind == 'pfree' else ([[y.min() * rng.uniform(0.05, 2.0)]] if self.kind == 'rcone' else [])
                p0 = np.concatenate([[K0], x0, Ld] + extra)
            lb = -np.inf * np.ones(self.npar()); ub = np.inf * np.ones(self.npar())
            if self.kind == 'pfree':
                lb[-1], ub[-1] = 0.3, 6.0
                p0[-1] = np.clip(p0[-1], 0.31, 5.9)
            if self.kind == 'rcone':
                lb[-1], ub[-1] = 0.0, 1.0
                p0[-1] = np.clip(p0[-1], 1e-6, 0.99)
            try:
                res = least_squares(lambda P: self.predict(P, X) - y, p0, bounds=(lb, ub), method='trf',
                                    x_scale='jac', max_nfev=20000)
            except Exception:
                continue
            if best is None or res.cost < best.cost:
                best = res
        return best
