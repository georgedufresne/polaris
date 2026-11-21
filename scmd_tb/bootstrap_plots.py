#!/usr/bin/env python3
"""Bootstrap figures: the null distributions, the star-count test, and realisations against their reference."""
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
OUT = f"{B}/bootstrap_plots"
H = "/dartfs-hpc/rc/home/8/f007058"
T = f"{B}/tuc_grid_nearest_bin02"
N6 = f"{B}/boot_ngc6791"
sys.path.insert(0, f"{T}/code_boot")
from isochronetoolbox.scmd import Setup, context                     # noqa: E402
from isochronetoolbox.scmd import bootstrap as BS                     # noqa: E402
from isochronetoolbox.scmd.pipeline import _window                    # noqa: E402
from isochronetoolbox.scmd.stats import ks_2d, rank2d                 # noqa: E402

INK, ACC, ACC2 = "#333333", "#2f6db5", "#c0392b"


def ctx_of(d):
    s = Setup(**{k: v for k, v in json.load(open(f"{d}/setup.json")).items() if k in Setup.__dataclass_fields__})
    return context(s), s


# ---------- 1. null distributions ------------------------------------------------------------
orig = pd.read_csv(f"{H}/Bootstrap/bootstrap_2d_ks_results.csv").ks_statistic.values
tb6 = pd.read_csv(f"{N6}/bootstrap_2d_ks_results.csv").ks_statistic.values
tuc = pd.read_csv(f"{T}/bootstrap_2d_ks_results.csv").ks_statistic.values
grid = pd.read_csv(f"{T}/combined_2d_ks_results.csv")
best_run = grid.groupby("simulation_id").ks_statistic.min().sort_values()

fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
bins = np.linspace(0, 0.075, 76)
ax[0].hist(orig, bins, density=True, histtype="step", lw=1.8, color=INK, label=f"original (Nov 2025): {orig.mean():.4f} +- {orig.std():.4f}")
ax[0].hist(tb6, bins, density=True, histtype="step", lw=1.8, color=ACC, label=f"toolbox: {tb6.mean():.4f} +- {tb6.std():.4f}")
for v, c in [(orig.mean() + 3 * orig.std(), INK), (tb6.mean() + 3 * tb6.std(), ACC)]:
    ax[0].axvline(v, color=c, ls="--", lw=1)
ax[0].set_title("NGC 6791 null: 10,000 realisations of run 20140, 8.4 Gyr\nvs its grid sCMD (dashed: mean + 3 sigma)", fontsize=10)
ax[0].set_xlabel("2D KS")
ax[0].set_ylabel("density")
ax[0].legend(fontsize=8)

b2 = np.linspace(0, 0.03, 121)
ax[1].hist(tuc, b2, density=True, color=ACC, alpha=0.85, label=f"null: {tuc.mean():.4f} +- {tuc.std():.4f}")
ax[1].axvline(tuc.mean() + 3 * tuc.std(), color=INK, ls="--", lw=1.2, label=f"mean + 3 sigma = {tuc.mean() + 3 * tuc.std():.4f}")
for k, v in enumerate(best_run.values[:10]):
    ax[1].axvline(v, color=ACC2, lw=1.0, alpha=0.9, label="best fit of the 10 best runs" if k == 0 else None)
ax[1].set_xlim(0, 0.03)
ax[1].set_title("47 Tuc null: 10,000 realisations of run 10099, 11.4 Gyr\nvs its grid sCMD, and the best fits against the data", fontsize=10)
ax[1].set_xlabel("2D KS")
ax[1].legend(fontsize=8)

# ---------- 2. the star-count test (47 Tuc subsampled to NGC 6791's 1,177 stars) --------------
ctx, s = ctx_of(T)
run, age, dm, e = 10099, 11.4, 13.38, 0.015
_, make = ctx.model(run, age)
sims = [make(np.random.RandomState(i), dm, e) for i in range(6)]
pool = np.vstack([ctx.points(c) for c in sims])
pool = pool[_window(pool, s.ks_window)]
obs = ctx.obs3
rs = np.random.RandomState(1)
N = 1177
def ks(a, b):
    return ks_2d(a, b, rank2d(a) / len(a))[0]
null_s, data_s = [], []
for _ in range(1000):
    null_s.append(ks(pool[rs.choice(len(pool), N, replace=False)], pool[rs.choice(len(pool), N, replace=False)]))
    data_s.append(ks(obs[rs.choice(len(obs), N, replace=False)], pool[rs.choice(len(pool), N, replace=False)]))
null_s, data_s = np.array(null_s), np.array(data_s)
lim = null_s.mean() + 3 * null_s.std()
b3 = np.linspace(0, 0.08, 81)
ax[2].hist(null_s, b3, density=True, histtype="step", lw=1.8, color=ACC, label=f"model vs model: {null_s.mean():.4f} +- {null_s.std():.4f}")
ax[2].hist(data_s, b3, density=True, histtype="step", lw=1.8, color=ACC2, label=f"data vs model: {data_s.mean():.4f} +- {data_s.std():.4f}")
ax[2].axvline(lim, color=INK, ls="--", lw=1.2, label=f"3 sigma {lim:.4f}: {np.mean(data_s <= lim) * 100:.0f}% of data draws pass")
ax[2].set_title(f"47 Tuc best cell with only {N} stars (NGC 6791's count):\nthe same misfit is lost in the noise", fontsize=10)
ax[2].set_xlabel("2D KS")
ax[2].legend(fontsize=8)
fig.tight_layout()
fig.savefig(f"{OUT}/bootstrap_distributions.png", dpi=130)
plt.close(fig)
print("star-count test:", null_s.mean(), null_s.std(), data_s.mean(), data_s.std(), np.mean(data_s <= lim))

# ---------- 3. realisations against the reference ------------------------------------------
def panel(ctx, s, cell, ids, path, labels, xlim, ylim, title):
    run, age, dm, e = cell
    ref = BS.reference(ctx, run, age, dm, e)
    ref_cdf = rank2d(ref) / len(ref)
    _, make = ctx.model(run, age)
    fig, axs = plt.subplots(1, len(ids) + 1, figsize=(3.6 * (len(ids) + 1), 5.2), sharex=True, sharey=True)
    small = len(ref) < 5000
    sz = 4 if small else 0.15
    axs[0].scatter(ref[:, 0], ref[:, 1], s=sz, c=INK, alpha=0.6 if small else 0.3, linewidths=0, rasterized=True)
    axs[0].set_title(f"reference: grid sCMD\n{len(ref)} stars in the KS window", fontsize=9)
    for a, i in zip(axs[1:], ids):
        cat = make(np.random.RandomState(i), dm, e)
        r = ctx.compare(ref, ref_cdf, cat)
        p = ctx.points(cat)
        p = p[_window(p, s.ks_window)]
        a.scatter(p[:, 0], p[:, 1], s=sz, c=ACC, alpha=0.6 if small else 0.3, linewidths=0, rasterized=True)
        a.set_title(f"bootstrap {i}\nKS {r['ks_statistic']:.4f}  JS {r['js_divergence']:.4f}", fontsize=9)
    for a in axs:
        a.set_xlim(*xlim)
        a.set_ylim(*ylim)
        a.set_xlabel(labels[0])
    axs[0].set_ylabel(labels[1])
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)

panel(ctx, s, (10099, 11.4, 13.38, 0.015), [0, 1, 2, 9999], f"{OUT}/tuc_realisations.png", ("F606W - F814W", "F606W"),
      (0.3, 1.0), (19.05, 14.95), "47 Tuc bootstrap: realisations of run 10099, 11.4 Gyr, mu 13.38, E(V-I) 0.015 against the grid sCMD")
c6, s6 = ctx_of(N6)
panel(c6, s6, (20140, 8.4, 13.33, 0.16), [0, 1, 2, 9999], f"{OUT}/ngc6791_realisations.png", ("BP - RP", "G"),
      (0.9, 1.7), (18.05, 16.2), "NGC 6791 bootstrap: realisations of run 20140, 8.4 Gyr, DM 13.33, E(B-V) 0.16 against the grid sCMD")
print("figures in", OUT)
