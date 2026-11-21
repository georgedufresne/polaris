#!/usr/bin/env python3
"""(M15 copy of best3_and_null.py) The three best runs of an M15 grid and where they sit in its bootstrap null.

    best3_and_null.py GRID_DIR [LABEL] [OUT.png]

Top row: the observed fit stars and the best sCMD of the three best runs (by 2D
KS, ages within each run's track coverage), regenerated with the grid's frozen
code.  Bottom row: the bootstrap null with the three fits marked, and the Hess
residuals (obs - sCMD) / sqrt(obs + sCMD).  Writes GRID_DIR/analysis/best3_and_null.png.
"""
import json, os, sys
G = sys.argv[1].rstrip("/")
LABEL = sys.argv[2] if len(sys.argv) > 2 else os.path.basename(G)
sys.path.insert(0, f"{G}/code")
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from isochronetoolbox.scmd import Setup, context, regenerate
from isochronetoolbox.scmd.pipeline import _window
B = "/dartfs-hpc/scratch/f007058/scmd_tb"
S = Setup(**{k: v for k, v in json.load(open(f"{G}/setup.json")).items() if k in Setup.__dataclass_fields__})
ctx = context(S)
lim_age = {int(k): v for k, v in json.load(open(os.environ.get("TUC_AGE_LIMITS", f"{B}/tuc_age_limits.json"))).items()}
g = pd.read_csv(f"{G}/combined_2d_ks_results.csv")
g = g[g.age_gyr < g.simulation_id.map(lim_age) - 1e-6]
best = g.loc[g.groupby("simulation_id").ks_statistic.idxmin()].nsmallest(3, "ks_statistic")
null = pd.read_csv(f"{G}/bootstrap_2d_ks_results.csv").ks_statistic.values
lim = null.mean() + 3 * null.std()
cell = json.load(open(f"{G}/bootstrap/cell.json"))


def mc(run):
    return dict((l.split()[1], float(l.split()[0])) for l in open(f"/dartfs-hpc/scratch/f007058/M15/dsep4/var/varfeh-227.{run}") if l.strip())


obs = ctx.obs3
zp = "zero-point corrected" if S.zp_correct else "no zero-point correction"
psf = "PSF scatter on" if (S.extra_sigma_table or S.extra_sigma_v or S.extra_sigma_i) else "no PSF scatter"
xb, yb = np.arange(0.32, 0.921, 0.02), np.arange(16.6, 20.701, 0.05)
ho, _, _ = np.histogram2d(obs[:, 0], obs[:, 1], [xb, yb])
C = ["#c0392b", "#e9a23b", "#6c5ce7"]
fig = plt.figure(figsize=(19, 10.5))
gs = GridSpec(2, 4, figure=fig, hspace=0.3, wspace=0.25)
ax0 = fig.add_subplot(gs[0, 0])
ax0.scatter(obs[:, 0], obs[:, 1], s=0.15, c="0.2", alpha=0.35, linewidths=0, rasterized=True)
ax0.set_title(f"Observed M15 ({zp})\n{len(obs)} fit stars", fontsize=10)
axes = [ax0]
for k, (r, c) in enumerate(zip(best.itertuples(), C), 1):
    cat = regenerate(ctx, int(r.simulation_id), r.age_gyr, r.distance_modulus, r.extinction_ebv)
    chk = ctx.stage3(cat)["ks_statistic"]
    p = ctx.points(cat); p = p[_window(p, S.ks_window)]
    m = mc(int(r.simulation_id))
    a = fig.add_subplot(gs[0, k], sharex=ax0, sharey=ax0); axes.append(a)
    a.scatter(p[:, 0], p[:, 1], s=0.15, c=c, alpha=0.35, linewidths=0, rasterized=True)
    a.set_title(f"#{k}  run {int(r.simulation_id)}: {r.age_gyr:.1f} Gyr, μ {r.distance_modulus:.2f}, E {r.extinction_ebv:.3f}\n"
                f"KS {r.ks_statistic:.4f}   [Fe/H] {m['FeH']:.2f}, Y {m['Yprim']:.3f}, α_MLT {m['CMIXLA']:.2f}", fontsize=9.5)
    hs, _, _ = np.histogram2d(p[:, 0], p[:, 1], [xb, yb])
    z = np.where(ho + hs > 0, (ho - hs) / np.sqrt(np.maximum(ho + hs, 1)), np.nan)
    ar = fig.add_subplot(gs[1, k], sharex=ax0, sharey=ax0); axes.append(ar)
    im = ar.pcolormesh(xb, yb, z.T, cmap="RdBu_r", vmin=-5, vmax=5, shading="flat")
    ar.set_title(f"#{k} residual (obs - sCMD) / sqrt(obs + sCMD)", fontsize=9.5)
    print(f"#{k} run {int(r.simulation_id)} {r.age_gyr:.1f} Gyr mu {r.distance_modulus:.2f} E {r.extinction_ebv:.3f} "
          f"KS {r.ks_statistic:.5f} (replay {chk:.5f}) p {(null >= r.ks_statistic).mean():.4f} "
          f"FeH {m['FeH']:.3f} Y {m['Yprim']:.3f} aMLT {m['CMIXLA']:.2f}")
for a in axes:
    a.set_xlim(0.36, 0.84); a.set_ylim(20.7, 16.55)
    a.set_xlabel("F606W - F814W")
ax0.set_ylabel("F606W")
fig.colorbar(im, ax=axes[2::2], shrink=0.75, pad=0.01, label="red: data has more stars   blue: sCMD has more")
an = fig.add_subplot(gs[1, 0])
hi = 1.12 * max(lim, best.ks_statistic.max(), np.percentile(null, 99.9))
an.hist(null, np.linspace(0, hi, 90), density=True, color="#2f6db5", alpha=0.85,
        label=f"bootstrap null, {len(null):,} draws\n{null.mean():.4f} ± {null.std():.4f}")
an.axvline(lim, color="#333333", ls="--", lw=1.5, label=f"mean + 3σ = {lim:.4f}")
for k, (r, c) in enumerate(zip(best.itertuples(), C), 1):
    an.axvline(r.ks_statistic, color=c, lw=2.2,
               label=f"#{k} run {int(r.simulation_id)}: {r.ks_statistic:.4f} (p = {(null >= r.ks_statistic).mean():.2f})")
an.set_xlim(0, hi)
an.set_xlabel("2D KS (95th percentile)"); an.set_ylabel("density")
an.set_title(f"bootstrap at run {cell['simulation_id']}, {cell['age_gyr']:.1f} Gyr: fresh stand-ins, resampled AST lists", fontsize=9)
an.legend(fontsize=7.5, loc="upper right")
fig.suptitle(f"M15 ({LABEL}: {zp}, {psf}): the three best fits and where they sit in the bootstrap null", fontsize=13)
out = sys.argv[3] if len(sys.argv) > 3 else f"{G}/analysis/best3_and_null.png"
os.makedirs(os.path.dirname(out), exist_ok=True)
fig.savefig(out, dpi=115, bbox_inches="tight")
print(out)
