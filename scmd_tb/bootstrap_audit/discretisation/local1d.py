"""Closed-form 3-point minima (hyperbola, parabola, symmetric V) across the trough, plus a 1D
hyperbola along mu at 11.4 Gyr / E(V-I)=0.015.  Values are single noisy draws (sigma ~0.0005-0.001)."""
import sys
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, '/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/discretisation')
from common import load


def three(km, k0, kp):
    """minima of y(t) through (-1, km), (0, k0), (+1, kp) for three shapes (t in steps)."""
    out = {}
    # hyperbola y^2 = K^2 + c^2 (t - t0)^2  (exact through 3 points)
    a = (km ** 2 + kp ** 2 - 2 * k0 ** 2) / 2          # c^2
    t0 = (km ** 2 - kp ** 2) / (4 * a) * -1 if a > 0 else np.nan
    t0 = -(kp ** 2 - km ** 2) / (4 * a)
    out['hyperbola'] = (np.sqrt(max(k0 ** 2 - a * t0 ** 2, 0)), t0)
    # parabola y = K + b (t - t0)^2
    b = (km + kp - 2 * k0) / 2
    t0p = -(kp - km) / (4 * b)
    out['parabola'] = (k0 - b * t0p ** 2, t0p)
    # symmetric V  y = K + s |t - t0|   (|t0| <= 1)
    B = (km - kp) / 2
    out['V'] = (k0 - abs(B), (B / ((km + kp) / 2 - (k0 - abs(B)))))
    return out


d = load()
g = d.groupby(['age', 'E']).ks.min()
cuts = {
    'age cut, E=0.020 (10.8/10.9/11.0)': (g[(10.8, 0.02)], g[(10.9, 0.02)], g[(11.0, 0.02)]),
    'age cut, E=0.015 (11.3/11.4/11.5)': (g[(11.3, 0.015)], g[(11.4, 0.015)], g[(11.5, 0.015)]),
    'age cut, E=0.010 (11.9/12.0/12.1)': (g[(11.9, 0.01)], g[(12.0, 0.01)], g[(12.1, 0.01)]),
    'E cut, 11.4 Gyr mu=13.38 (0.010/0.015/0.020)': tuple(d[(d.src == 'allcells') & np.isclose(d.mu, 13.38)].set_index('E').ks.loc[[0.01, 0.015, 0.02]]),
    'E cut, 11.4 Gyr best mu per E (0.010/0.015/0.020)': (g[(11.4, 0.01)], g[(11.4, 0.015)], g[(11.4, 0.02)]),
}
for name, (km, k0, kp) in cuts.items():
    r = three(km, k0, kp)
    print(f'{name}: values {km:.4f} {k0:.4f} {kp:.4f}')
    for k, (K, t0) in r.items():
        print(f'    {k:9s} min {K:.5f} at {t0:+.3f} steps   (node minus min = {k0 - K:.5f})')
    # same with the node replaced by its expected value 0.0178 (fresh realisations) for the 11.4 node
    if '11.3/11.4/11.5' in name or 'mu=13.38' in name:
        r2 = three(km, 0.0178, kp)
        print('    with the node at its fresh-realisation mean 0.0178: ' +
              ', '.join(f'{k} {K:.5f}' for k, (K, t0) in r2.items()))

# 1D hyperbola along mu at 11.4 Gyr, E = 0.015
a = d[(d.src == 'allcells') & np.isclose(d.E, 0.015)].sort_values('mu')
for T in (0.03, 0.04, 0.08):
    s = a[a.ks < T]
    x, y = s.mu.to_numpy(), s.ks.to_numpy()
    f = lambda p: np.sqrt(p[0] ** 2 + p[1] ** 2 * (x - p[2]) ** 2) - y
    r = least_squares(f, [0.016, 0.25, 13.39])
    K, c, m0 = abs(r.x[0]), abs(r.x[1]), r.x[2]
    nodes = np.round(np.arange(13.15, 13.455, 0.01), 2)
    pn = np.sqrt(K ** 2 + c ** 2 * (nodes - m0) ** 2)
    print(f'mu profile (KS<{T}, {len(s)} cells): K_mu={K:.5f}  mu*={m0:.4f}  slope={c:.3f}/mag  rms={np.sqrt(np.mean(r.fun**2)):.5f}  '
          f'best-node pred={pn.min():.5f} (mu-lattice penalty {pn.min() - K:.6f}); pred at 13.38: {np.sqrt(K**2 + c**2*(13.38-m0)**2):.5f}')
