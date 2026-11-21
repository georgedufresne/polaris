#!/usr/bin/env python3
"""Figures for the bootstrap / artificial-star error model work (47 Tuc)."""
import json, sys
from dataclasses import replace
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
OUT = f"{B}/bootstrap_fix_plots"
sys.path.insert(0, f"{B}/smooth_dev")
from isochronetoolbox.scmd import Setup, context                      # noqa: E402
from isochronetoolbox.scmd.ast import NearestErrors, ArtificialStars, centre  # noqa: E402

C = ["#2f6db5", "#c0392b", "#2a9d8f", "#e9a23b", "#6c5ce7"]   # fixed categorical order
INK, GREY = "#333333", "#999999"
plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True,
                     "grid.alpha": 0.25, "font.size": 10})

G = f"{B}/tuc_grid_nearest_bin02"
S0 = Setup(**{k: v for k, v in json.load(open(G + "/setup.json")).items() if k in Setup.__dataclass_fields__})

# ------------------------------------------------------------- figure 1: the nulls
nulls = [("copied errors, one fixed list (first run)", f"{G}/bootstrap_2d_ks_results.csv", "ks_statistic", C[0]),
         ("copied errors, independent halves (removed)", f"{B}/tuc_grid_nearest_bin02_split/bootstrap_2d_ks_results.csv", "ks_statistic", GREY),
         ("smoothed errors, one fixed list", f"{B}/smooth_test/smooth_noresample/bootstrap_2d_ks_results.csv", "ks_statistic", C[2]),
         ("smoothed errors (100 nbrs), resampled list per realisation", f"{B}/smooth_test/smooth_resample/bootstrap_2d_ks_results.csv", "ks_statistic", C[1]),
         ("NEW MODEL, one stand-in, resampled lists", f"{B}/smooth_test/newmodel/bootstrap_2d_ks_results.csv", "ks_statistic", C[3]),
         ("NEW MODEL, fresh stand-in every draw, resampled lists", f"{B}/smooth_test/newmodel_fresh/bootstrap_2d_ks_results.csv", "ks_statistic", C[4])]
dvm = pd.read_csv(f"{B}/smooth_test/smooth_resample/data_vs_model.csv").ks_statistic
dvn = pd.read_csv(f"{B}/smooth_test/newmodel/data_vs_model.csv").ks_statistic
fig, ax = plt.subplots(1, 2, figsize=(15, 5.2))
bins = np.linspace(0, 0.022, 89)
for lab, path, col, c in nulls:
    v = pd.read_csv(path)[col].values
    ax[0].hist(v, bins, density=True, histtype="step", lw=2, color=c,
               label=f"{lab}: {v.mean():.4f} ± {v.std():.4f}, 3σ {v.mean() + 3 * v.std():.4f}")
ax[0].axvspan(dvm.mean() - dvm.std(), dvm.mean() + dvm.std(), color=INK, alpha=0.15, lw=0)
ax[0].axvline(dvm.mean(), color=INK, lw=1.5, label=f"data vs model (smoothed, 100 nbrs), 12 realisations: {dvm.mean():.4f} ± {dvm.std():.4f}")
ax[0].axvspan(dvn.mean() - dvn.std(), dvn.mean() + dvn.std(), color=C[4], alpha=0.15, lw=0)
ax[0].axvline(dvn.mean(), color=C[4], lw=1.8, label=f"data vs NEW MODEL, 12 realisations: {dvn.mean():.4f} ± {dvn.std():.4f}")
ax[0].axvline(0.0168, color=INK, lw=1.2, ls="--", label="best grid fit (copied errors): 0.0168")
ax[0].set_xlabel("2D KS (95th percentile)")
ax[0].set_ylabel("density")
ax[0].set_title("47 Tuc bootstrap nulls at the best cell (run 10099, 11.4 Gyr, μ 13.38, E(V-I) 0.015)", fontsize=10)
ax[0].legend(fontsize=7.5, loc="upper right")
y = pd.read_csv(f"{B}/smooth_test/ying_null/ying_null.csv")
bins2 = np.linspace(0, 0.07, 141)
ref = pd.read_csv(nulls[5][1]).ks_statistic.values
ax[1].hist(ref, bins2, density=True, histtype="step", lw=2, color=C[4], label=f"ours: new model, fresh stand-in, resampled lists: {ref.mean():.4f}")
for (cfg, size), c in zip([("binned per-band (M. Ying)", "10 x N_obs"), ("smoothed joint (ours)", "10 x N_obs"),
                           ("smoothed joint (ours)", "1 x N_obs")], [C[4], C[3], C[2]]):
    g = y[(y.config == cfg) & (y["size"] == size)].ks
    ax[1].hist(g, bins2, density=True, histtype="stepfilled", alpha=0.55, color=c,
               label=f"M. Ying's resample bootstrap, {cfg.replace(' (M. Ying)', '').replace(' (ours)', '')} errors, {size}: {g.mean():.4f}")
ax[1].axvline(0.0168, color=INK, lw=1.2, ls="--", label="best grid fit: 0.0168")
ax[1].set_xlabel("2D KS (95th percentile)")
ax[1].set_title("M. Ying's resample bootstrap (observed stars resampled, artificial-star errors added again)\n"
                "measures the extra broadening of a second error layer, not sampling noise", fontsize=10)
ax[1].legend(fontsize=7.5, loc="upper right")
fig.tight_layout()
fig.savefig(f"{OUT}/1_bootstrap_nulls.png", dpi=130)
plt.close(fig)

# ------------------------------------------------------------- figure 2: the artificial-star error model
a = pd.read_csv(S0.ast, dtype=np.float32)
vin, iin = a.InputF606W + 30.5742, a.InputF814W + 29.6777
dv, di = a.InputF606W - a.OutputF606W, a.InputF814W - a.OutputF814W
rec = (np.hypot(a.InputX - a.OutputX, a.InputY - a.OutputY) < 0.5) & (dv.abs() < 0.75) & (di.abs() < 0.75)
xc, yc = centre(a.InputX.values, a.InputY.values)
r = np.hypot(a.InputX - xc, a.InputY - yc)
mad = lambda x: 1.4826 * np.median(np.abs(x - np.median(x)))
fig, ax = plt.subplots(2, 3, figsize=(16, 9.5))
redges = np.array([0, 150, 300, 450, 600, 800, 1000, 1300, 1700, 2200, 3000])
rc_ = 0.5 * (redges[1:] + redges[:-1])
for (lo, hi), c in zip([(16.8, 17.3), (18.3, 18.8)], [C[0], C[1]]):
    s = (vin > lo) & (vin < hi)
    sig = [mad(dv[s & rec & (r >= a0) & (r < a1)]) for a0, a1 in zip(redges[:-1], redges[1:])]
    comp = [rec[s & (r >= a0) & (r < a1)].mean() for a0, a1 in zip(redges[:-1], redges[1:])]
    ax[0, 0].plot(rc_, sig, "o-", color=c, lw=2, label=f"F606W {lo}-{hi}")
    ax[0, 1].plot(rc_, comp, "o-", color=c, lw=2, label=f"F606W {lo}-{hi}")
ax[0, 0].set_yscale("log"); ax[0, 0].set_xlabel("distance from the centre (px)"); ax[0, 0].set_ylabel("error scale, 1.48 MAD of ΔV (mag)")
ax[0, 0].set_title("errors fall 30-40% per 150 px inside 800 px,\n~20% per 350 px outside", fontsize=10); ax[0, 0].legend(fontsize=8)
ax[0, 1].set_xlabel("distance from the centre (px)"); ax[0, 1].set_ylabel("completeness"); ax[0, 1].legend(fontsize=8)
ax[0, 1].set_title("completeness: faint stars lose up to 26% in the core", fontsize=10)
s = (r > 800) & (r < 1300) & rec
vb = np.arange(15.0, 19.25, 0.25)
sv = [mad(dv[s & (vin >= v0) & (vin < v0 + 0.25)]) for v0 in vb]
ax[0, 2].plot(vb + 0.125, sv, "o-", color=C[2], lw=2)
ax[0, 2].axvline(16.48, color=INK, ls="--", lw=1)
ax[0, 2].text(16.52, max(sv) * 0.95, "short | long exposures\n(F606W 16.48)", fontsize=8, va="top")
ax[0, 2].set_xlabel("input F606W"); ax[0, 2].set_ylabel("error scale (mag), radius 800-1300 px")
ax[0, 2].set_title("errors double across the exposure boundary:\nneighbourhoods must not straddle it", fontsize=10)
# copied vs smoothed draws at a crowded position, and the neighbour spread of the two metrics
copy = NearestErrors(a, k=30)
smooth_lin = NearestErrors(a, k=30, smooth=True)
split = (16.48 - 30.5742, 15.55 - 29.6777)
crowd = NearestErrors(a, k=30, metric="crowding", split=split)
m = 60000
args = (np.full(m, 300.0), np.full(m, 18.5 - 30.5742), np.full(m, 18.5 - 0.58 - 29.6777))
oc, evc, eic = copy.draw(np.random.RandomState(1), *args)
os_, evs, eis = smooth_lin.draw(np.random.RandomState(1), *args)
bb = np.linspace(-0.25, 0.25, 101)
ax[1, 0].hist(evc[oc], bb, density=True, histtype="step", lw=2, color=C[0], label=f"copied ({len(np.unique(evc[oc]))} distinct values)")
ax[1, 0].hist(evs[os_], bb, density=True, histtype="step", lw=2, color=C[1], label=f"smoothed ({len(np.unique(evs[os_]))} distinct values)")
ax[1, 0].set_yscale("log"); ax[1, 0].set_xlabel("ΔV = input - output (mag)"); ax[1, 0].set_ylabel("density")
ax[1, 0].set_title("errors drawn for a star at F606W 18.5, r = 300 px", fontsize=10); ax[1, 0].legend(fontsize=8)
ax[1, 1].hist((evc - eic)[oc], bb / 2, density=True, histtype="step", lw=2, color=C[0], label=f"copied: corr(ΔV, ΔI) {np.corrcoef(evc[oc], eic[oc])[0, 1]:.2f}")
ax[1, 1].hist((evs - eis)[os_], bb / 2, density=True, histtype="step", lw=2, color=C[1], label=f"smoothed: corr {np.corrcoef(evs[os_], eis[os_])[0, 1]:.2f}")
ax[1, 1].set_yscale("log"); ax[1, 1].set_xlabel("colour error Δ(V-I) (mag)")
ax[1, 1].set_title("smoothing keeps the colour error and the band correlation", fontsize=10); ax[1, 1].legend(fontsize=8)
for ne, lab, c in [(copy, "linear metric (0.1 mag, 200 px)", C[0]), (crowd, "crowding metric (0.1 mag, ln(r+200)/0.09)", C[1])]:
    for R, ls in [(250.0, "-"), (1500.0, ":")]:
        q = np.random.RandomState(5)
        vv = q.uniform(17.0, 18.8, 3000)
        pts = ne._coords(vv - 30.5742, vv - 0.57 - 29.6777, np.full(3000, R))
        _, near = ne.tree.query(pts)
        rr = np.hypot(a.InputX.values - xc, a.InputY.values - yc)[ne.nbr[near]]
        spread = np.percentile(np.abs(rr - R), 84, axis=1)
        ax[1, 2].hist(spread, np.linspace(0, 400, 81), histtype="step", lw=2, ls=ls, color=c,
                      label=f"{lab}, star at r = {int(R)} px")
ax[1, 2].set_xlabel("84th percentile of |Δr| over the 30 neighbours (px)"); ax[1, 2].set_ylabel("stars")
ax[1, 2].set_title("how far in radius the chosen neighbours reach", fontsize=10); ax[1, 2].legend(fontsize=7.5)
fig.suptitle("47 Tuc artificial-star error model: what the errors depend on, and copying vs smoothed sampling", fontsize=12)
fig.tight_layout()
fig.savefig(f"{OUT}/2_ast_error_model.png", dpi=125)
plt.close(fig)

# ------------------------------------------------------------- figure 3: completeness-corrected radii
ctx = context(replace(S0, ast_errors="smooth", ast_k=30, ast_metric="crowding", radial="corrected", crn=True))
A = ctx.ast
fig, ax = plt.subplots(1, 2, figsize=(13, 4.8))
for vb_, c in zip([15.6, 17.1, 18.6], [C[0], C[2], C[1]]):
    k = np.clip(np.digitize(vb_, A._redges) - 1, 0, A._rbins - 1)
    ax[0].plot(A._rdist[k], A._rcdf[k], color=c, lw=2, label=f"corrected, F606W {A._redges[k]:.2f}-{A._redges[k + 1]:.2f} (mean completeness {A._rcomp[k]:.3f})")
ax[0].plot(A.dist, A.cdf, color=INK, lw=1.5, ls="--", label="observed, all magnitudes (what the grid used)")
ax[0].set_xlim(0, 1500); ax[0].set_xlabel("distance from the centre (px)"); ax[0].set_ylabel("cumulative fraction")
ax[0].set_title("radial distributions the sCMD stars are drawn from", fontsize=10); ax[0].legend(fontsize=7.5)
frac = [np.interp(400, A._rdist[k], A._rcdf[k]) for k in range(A._rbins)]
cen = 0.5 * (A._redges[1:] + A._redges[:-1])
ax[1].plot(cen, frac, "o-", color=C[1], lw=2, label="corrected, per magnitude bin")
ax[1].axhline(np.interp(400, A.dist, A.cdf), color=INK, ls="--", lw=1.5, label="observed, all magnitudes")
ax[1].set_xlabel("F606W bin"); ax[1].set_ylabel("fraction of stars within 400 px")
ax[1].set_title("core concentration by magnitude: incompleteness undone,\nmass segregation kept", fontsize=10); ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(f"{OUT}/3_radial_correction.png", dpi=130)
plt.close(fig)
print("figures in", OUT)
