"""(c) How far is the grid minimum (min over 97 runs x 61 ages x 4 cells) biased low?

y_c = mu_c + sigma_c z_c: grid KS of cell c = its expected KS (over model realisations,
data fixed) + realisation noise.  sigma_c measured from fresh realisations of five top cells.

1. direct: grid value vs fresh-realisation mean, for the five measured cells
2. empirical Bayes (g-modelling): deconvolve the distribution g(mu) of expected KS from the
   grid's y (cells with y < YMAX) -- NPMLE (EM on a fine grid) and a smooth log-polynomial
   g -- then
   * posterior mean E[mu | y] for the measured cells (validation against the fresh means;
     selection on y does not bias a Bayes posterior when g is right)
   * posterior-predictive simulation of the whole grid: mu_c ~ p(mu | y_c), y*_c = mu_c +
     sigma_c z: E[mu(argmin) - min y*] (winner's curse of the selected cell) and
     E[min mu - min y*] (minimum vs the truly best expected KS)
3. plug-in check: y*_c = y_c + sigma z (no deconvolution)
Sensitivity: sigma = 0.0005, 0.0007, 0.0008, 0.00095, and a heteroscedastic mixture of the
five measured sigmas.
"""
import numpy as np
import pandas as pd
from scipy.optimize import minimize

O = "/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/reference-selection/out"
grid = pd.read_csv("/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02/combined_2d_ks_results.csv")
y_all = grid.ks_statistic.values
YMAX = 0.035
y = y_all[y_all < YMAX]
print(f"grid cells: {len(y_all)}; used (y < {YMAX}): {len(y)}; min {y.min():.6f}")

# ---- measured cells -------------------------------------------------------------
sh = pd.read_csv(f"{O}/shared_design_all.csv")
oc = pd.read_csv(f"{O}/other_cells.csv")
cells = [("10099 11.4 13.38 0.015 (selected)", 0.0168281775635039, sh.data_ks.values)]
for (run, age, mu, e), g in oc.groupby(["run", "age", "mu", "evi"], sort=False):
    cells.append((f"{run} {age} {mu} {e}", g.grid_ks.iloc[0], g.data_ks.values))
print("\n1. DIRECT: grid KS vs data-vs-fresh-realisation KS (latest code, full AST list)")
sig_meas = []
for name, gk, x in cells:
    m, s, n = x.mean(), x.std(ddof=1), len(x)
    sig_meas.append(s)
    print(f"  {name:36s} grid {gk:.6f}  fresh {m:.6f} +- {s / np.sqrt(n):.6f} (sd {s:.6f}, n {n})  "
          f"grid - fresh {gk - m:+.6f} ({(gk - m) / s:+.2f} sd)")
sig_meas = np.array(sig_meas)
dfs = np.array([len(x) - 1 for _, _, x in cells])
sig_pool = np.sqrt(np.sum(dfs * sig_meas ** 2) / dfs.sum())
from scipy.stats import bartlett
print(f"  pooled sigma {sig_pool:.6f}; Bartlett test of equal sigmas p = {bartlett(*[x for _, _, x in cells]).pvalue:.3g}")
print(f"  mean (grid - fresh) over the 5 cells {np.mean([gk - x.mean() for _, gk, x in cells]):+.6f}")

# ---- g-modelling ------------------------------------------------------------------
support = np.round(np.arange(0.0130, YMAX + 0.0010 + 1e-9, 0.0001), 6)


def lik(yv, sig):
    """L[c, i] = p(y_c | mu_i); sig scalar or (values, probs) mixture."""
    if np.isscalar(sig):
        return np.exp(-0.5 * ((yv[:, None] - support[None, :]) / sig) ** 2) / sig
    vals, probs = sig
    return sum(p * np.exp(-0.5 * ((yv[:, None] - support[None, :]) / s) ** 2) / s for s, p in zip(vals, probs))


def trunc_norm(sig):
    """P(y < YMAX | mu_i): the used cells are truncated at YMAX."""
    from scipy.stats import norm
    if np.isscalar(sig):
        return norm.cdf((YMAX - support) / sig)
    vals, probs = sig
    return sum(p * norm.cdf((YMAX - support) / s) for s, p in zip(vals, probs))


def npmle(L, T, iters=3000):
    w = np.full(L.shape[1], 1.0 / L.shape[1])
    for _ in range(iters):
        # truncated EM: weights of the observed (y < YMAX) population
        f = L @ (w * 1.0)
        w = w * (L.T @ (1.0 / f)) / L.shape[0]
        w = w / w.sum()
    # w is the distribution of mu among cells that were observed below YMAX, as the
    # likelihood L already conditions on nothing; T is kept for reference
    return w


def logpoly(L, deg=6, lam=1e-3):
    xs = (support - support.mean()) / (support.max() - support.min())
    Q = np.vstack([xs ** k for k in range(1, deg + 1)]).T

    def nll(b):
        a = Q @ b
        w = np.exp(a - a.max())
        w /= w.sum()
        return -np.sum(np.log(L @ w + 1e-300)) + lam * np.sum(b ** 2)
    r = minimize(nll, np.zeros(deg), method="BFGS", options={"maxiter": 2000})
    a = Q @ r.x
    w = np.exp(a - a.max())
    return w / w.sum()


rng = np.random.default_rng(12345)
SIGMAS = [0.0005, 0.0007, 0.0008, 0.00095, ("mix", (sig_meas, np.full(len(sig_meas), 1 / len(sig_meas))))]
measured_y = np.array([gk for _, gk, _ in cells])
measured_m = np.array([x.mean() for _, _, x in cells])
measured_se = np.array([x.std(ddof=1) / np.sqrt(len(x)) for _, _, x in cells])
rows = []
print("\n2. EMPIRICAL BAYES (g-modelling) + posterior-predictive simulation of the grid minimum")
for sig in SIGMAS:
    label = f"sigma {sig:.5f}" if np.isscalar(sig) else "sigma mix of the 5 measured"
    s_ = sig if np.isscalar(sig) else sig[1]
    L = lik(y, s_)
    for gname, w in [("NPMLE", npmle(L, None)), ("logpoly6", logpoly(L))]:
        # posterior means of the measured cells
        Lm = lik(measured_y, s_)
        post = Lm * w[None, :]
        post /= post.sum(axis=1, keepdims=True)
        pm = post @ support
        # posterior-predictive simulation of the grid
        P = L * w[None, :]
        P /= P.sum(axis=1, keepdims=True)
        C = np.cumsum(P, axis=1)
        nsim = 2000
        wc, bm, same, ymin_s = [], [], [], []
        for _ in range(nsim):
            u = rng.random(len(y))
            idx = np.minimum((C < u[:, None]).sum(axis=1), len(support) - 1)
            mu = support[idx] + rng.uniform(-5e-5, 5e-5, len(y))      # jitter within the support step
            if np.isscalar(s_):
                sc = np.full(len(y), s_)
            else:
                sc = rng.choice(s_[0], len(y), p=s_[1])
            ys = mu + sc * rng.standard_normal(len(y))
            k = np.argmin(ys)
            wc.append(mu[k] - ys[k])
            bm.append(mu.min() - ys[k])
            same.append(k == np.argmin(mu))
            ymin_s.append(ys[k])
        wc, bm = np.array(wc), np.array(bm)
        rows.append({"sigma": label, "g": gname, "E[mu_sel - min y]": wc.mean(), "sd": wc.std(),
                     "E[min mu - min y]": bm.mean(), "sd2": bm.std(), "P(selected = best)": np.mean(same),
                     "E[min y*]": np.mean(ymin_s),
                     "post mean selected": pm[0], **{f"post mean cell{j + 1}": pm[j] for j in range(1, len(pm))},
                     "chi2 vs measured (5 cells)": float(np.sum(((pm - measured_m) / measured_se) ** 2))})
        print(f"  {label:28s} {gname:8s}: winner's curse E[mu_sel - min y] = {wc.mean():.6f} (sd {wc.std():.6f}); "
              f"E[min mu - min y] = {bm.mean():+.6f} (sd {bm.std():.6f}); P(selected cell is best) {np.mean(same):.2f}; "
              f"E[min y*] {np.mean(ymin_s):.6f}")
        print(f"      posterior means of measured cells: " + " ".join(f"{v:.5f}" for v in pm)
              + f"   (measured: " + " ".join(f"{v:.5f}" for v in measured_m) + ")")
pd.DataFrame(rows).to_csv(f"{O}/analyse_c_eb.csv", index=False)

# ---- plug-in check -------------------------------------------------------------------
print("\n3. PLUG-IN (mu = observed y, no deconvolution; spreads the near-ties, so a lower bound)")
for sig in [0.0005, 0.0007, 0.0008, 0.00095]:
    d = []
    for _ in range(4000):
        ys = y + sig * rng.standard_normal(len(y))
        d.append(y.min() - ys.min())
    print(f"  sigma {sig:.5f}: E[min y - min(y + noise)] = {np.mean(d):.6f} (sd {np.std(d):.6f})")

# ---- order-statistics intuition -------------------------------------------------------
print("\n4. ORDER STATISTICS: K cells tied in expectation, sigma 0.0008-0.00095: min biased by sigma * E[max of K normals]")
z = rng.standard_normal((20000, 64))
for K in (1, 2, 4, 8, 16, 32, 64):
    e = np.mean(z[:, :K].max(axis=1))
    print(f"  K {K:3d}: E[max] {e:.3f}  -> bias {0.0008 * e:.6f} .. {0.00095 * e:.6f}")
