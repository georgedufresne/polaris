"""Tables only: KS surface around the best cell (off-grid penalty), selection, 0.0468 filter."""
import numpy as np, pandas as pd
G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
SP = G + '_split'
c = pd.read_csv(G + '/all_cells_10099_11.4.csv')
print(c.columns.tolist(), len(c))
b = c.loc[c.ks.idxmin()]; print("best of 341 cells:", b.to_dict())
# quadratic fit of KS around the minimum: KS = a + q(dmu, de) using cells within +-0.04 mu, +-0.02 E
near = c[(np.abs(c.mu - 13.38) <= 0.0401) & (np.abs(c.evi - 0.015) <= 0.0201)].copy()
x = near.mu.values - 13.38; y = near.evi.values - 0.015
A = np.column_stack([np.ones_like(x), x, y, x * x, x * y, y * y])
coef, *_ = np.linalg.lstsq(A, near.ks.values, rcond=None)
res = near.ks.values - A @ coef
print("quadratic fit on", len(near), "cells; rms residual %.5f" % res.std())
a0, ax, ay, axx, axy, ayy = coef
H = np.array([[2 * axx, axy], [axy, 2 * ayy]])
xm = np.linalg.solve(H, -np.array([ax, ay]))
kmin = a0 + 0.5 * np.array([ax, ay]) @ xm
print("fitted minimum at dmu %.4f dE %.4f, KS_min %.5f; Hessian diag %.2f %.2f off %.2f" % (xm[0], xm[1], kmin, H[0, 0], H[1, 1], H[0, 1]))
# expected excess for a true position uniform within +-half a grid step around the nearest node
rng = np.random.RandomState(0)
u = rng.uniform(-0.5, 0.5, (200000, 2)) * np.array([0.01, 0.005])
exc = 0.5 * np.einsum('ni,ij,nj->n', u, H, u)
print("off-grid excess in KS (quadratic, uniform within half-step box): mean %.5f  max %.5f" % (exc.mean(), exc.max()))
# but KS is not quadratic near the minimum when noise-floor dominated: add in quadrature form too
# KS(dmu) ~ sqrt(KS0^2 + (g*dmu)^2): estimate slopes from the neighbours along each axis
for ax_name, sub in (("mu", c[np.isclose(c.evi, 0.015)].sort_values("mu")), ("evi", c[np.isclose(c.mu, 13.38)].sort_values("evi"))):
    print(ax_name, "slice:", " ".join(f"{v:.3f}:{k:.4f}" for v, k in zip(sub[ax_name if ax_name == 'mu' else 'evi'].values, sub.ks.values) if abs(v - (13.38 if ax_name == 'mu' else 0.015)) < 0.051))

comb = pd.read_csv(G + '/combined_2d_ks_results.csv')
print("stage-3 rows", len(comb), " runs", comb.simulation_id.nunique(), " ages", comb.age_gyr.nunique())
print("min KS row:", comb.loc[comb.ks_statistic.idxmin(), ["simulation_id", "age_gyr", "distance_modulus", "extinction_ebv", "ks_statistic", "kde_rank"]].to_dict())
print("rows with KS <= 0.0468 (setup ks_3sigma, NGC 6791 value):", int((comb.ks_statistic <= 0.0468).sum()), "of", len(comb))
f = pd.read_csv(G + '/filtered_ks_scmd.csv'); print("filtered_ks_scmd.csv rows:", len(f))
for lim, nm in ((0.007434816597419277, "shared 3sig"), (0.014632839607227669, "split 3sig")):
    print(f"rows with KS <= {lim:.4f} ({nm}):", int((comb.ks_statistic <= lim).sum()))
# best per (run, age), how many below various levels; distribution of the lowest 20
best = comb.sort_values("ks_statistic").head(20)
print(best[["simulation_id", "age_gyr", "distance_modulus", "extinction_ebv", "ks_statistic", "kde_rank"]].to_string(index=False))
# the null files
for nm, path in (("shared", G + '/bootstrap_2d_ks_results.csv'), ("split", SP + '/bootstrap_2d_ks_results.csv')):
    z = pd.read_csv(path)
    k = z.ks_statistic.values
    print(f"{nm}: n {len(k)} mean {k.mean():.5f} sd {k.std():.5f} p95 {np.percentile(k, 95):.5f} p99.87 {np.percentile(k, 99.865):.5f} max {k.max():.5f}"
          f"  n_obs uniq {z.n_observed.unique()[:3]} n_sim uniq {z.n_simulated.unique()[:3]}  ids 0..{z.bootstrap_number.max()} dup {z.bootstrap_number.duplicated().sum()}")
s1 = pd.read_csv(G + '/bootstrap_2d_ks_results.csv').set_index("bootstrap_number").ks_statistic
s2 = pd.read_csv(SP + '/bootstrap_2d_ks_results.csv').set_index("bootstrap_number").ks_statistic
print("corr(shared_i, split_i) over the same seeds:", float(np.corrcoef(s1.loc[s2.index], s2)[0, 1]))
