"""W14feh copy of ../scmd_W14_pipeline/age_breakdown.py: the W v1.4 leg is 2d_ks_W14feh (W fed [Fe/H] for our Z
convention); home inputs read from inputs/ (md5-checked copies). Writes age_breakdown.txt/png in this folder.
"""
#!/usr/bin/env python3
"""NGC 6791 sCMD age breakdown: W v1.4 (2d_ks_W14, interim while the grid runs) vs OG vs DSEP3 with the physics on
("New", analysis_noco). Same KS cut (0.0468) and weights (KDE of the original bootstrap null at each KS value, one best
DM/E per (MC set, age)) as scmd_3way_compare.py / scmd_5way_compare.py.

1. Posterior age per leg: weighted mean, sd, 2.5/16/50/84/97.5% quantiles, mode; on the MC sets all three keep, and on each
   leg's own sets.
2. Run by run (sets both legs keep): the KS-weighted mean age of each MC set, and its best-fit (lowest-KS) age; the
   W v1.4 - DSEP3 difference per set (median, scatter, fraction within 0.25 / 0.5 Gyr), and New - OG for reference.
3. Which inputs drive W v1.4 - New per set: Spearman correlation with each var-file input (initial_vars.csv).
4. Fit quality per set: lowest KS, W v1.4 vs New.
5. By surface boundary condition (KTTAU 0 / 1 / 5): per-set differences and posteriors on the common sets.
Writes age_breakdown.txt and age_breakdown.png here."""
import glob
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt     # noqa: E402
import numpy as np                  # noqa: E402
import pandas as pd                 # noqa: E402
from scipy.stats import gaussian_kde, spearmanr   # noqa: E402

S, H = "/dartfs-hpc/scratch/f007058", "/dartfs-hpc/rc/home/8/f007058"
HERE = os.path.dirname(os.path.abspath(__file__))
CUT = 0.0468
boot = pd.read_csv(f"{S}/scmd_W14feh_pipeline/inputs/bootstrap_2d_ks_results.csv").ks_statistic.values
kde = gaussian_kde(boot, bw_method="scott")
bins = np.arange(3.95, 12.15, 0.1)
ctr = (bins[:-1] + bins[1:]) / 2
iv = pd.read_csv(f"{S}/analysis_W/initial_vars.csv").rename(columns={"file_id": "simulation_id"}).set_index("simulation_id")


KEEP = set(iv.index[iv.KTTAU.isin([0, 5])])     # W v1.4 is not trained on KTTAU = 1: every leg restricted to KTTAU 0/5


def prep(d):
    d = d[d.simulation_id.isin(KEEP)]
    raw_best = d.loc[d.groupby("simulation_id").ks_statistic.idxmin(), ["simulation_id", "ks_statistic"]].set_index("simulation_id").ks_statistic
    d = d[d.ks_statistic <= CUT].copy()
    d["w"] = kde(d.ks_statistic.values)
    d = d.loc[d.groupby(["simulation_id", "age_gyr"]).ks_statistic.idxmin()].reset_index(drop=True)
    return d, raw_best


w14raw = pd.concat([pd.read_csv(f) for f in glob.glob(f"{S}/2d_ks_W14feh/2d_ks_results_*.csv")], ignore_index=True)
L = {"OG": prep(pd.read_csv(f"{S}/scmd_W14feh_pipeline/inputs/filtered_ks_scmd_OG.csv")),
     "DSEP3 physics on": prep(pd.read_csv(f"{S}/analysis_noco/filtered_ks_scmd.csv")),
     "W v1.4 Z-matched": prep(w14raw)}
NAMES = list(L)
COL = {"OG": "#52514e", "DSEP3 physics on": "#2a78d6", "W v1.4 Z-matched": "#eb6834"}
out = []
A = out.append


def post(d):
    h, _ = np.histogram(d.age_gyr, bins=bins, weights=d.w)
    p = h / h.sum()
    m = (ctr * p).sum()
    sd = np.sqrt(((ctr - m) ** 2 * p).sum())
    c = np.cumsum(p)
    q = [np.interp(x, c, ctr) for x in (0.025, 0.16, 0.5, 0.84, 0.975)]
    return m, sd, q, ctr[np.argmax(p)], p


def per_set(d):
    g = d.groupby("simulation_id")
    wm = pd.Series({k: np.average(x.age_gyr, weights=x.w) for k, x in g})
    best = d.loc[g.ks_statistic.idxmin()].set_index("simulation_id").age_gyr
    return pd.DataFrame({"wmean": wm, "best": best, "wsum": g.w.sum()})


n_w14 = w14raw.simulation_id.nunique()
A(f"NGC 6791 sCMD age breakdown -- W v1.4 with [Fe/H] for our Z convention ({n_w14} runs with 2D-KS results), OG, DSEP3 physics on")
A("=" * 112)
common = set.intersection(*[set(v[0].simulation_id) for v in L.values()])
for title, sel in (("same MC sets (kept by all three)", {k: v[0][v[0].simulation_id.isin(common)] for k, v in L.items()}),
                   ("each leg's own sets", {k: v[0] for k, v in L.items()})):
    A(f"\n1. Posterior age, {title}")
    A(f"{'leg':<18}{'sets':>6}{'mean':>8}{'sd':>7}{'2.5%':>8}{'16%':>8}{'median':>8}{'84%':>8}{'97.5%':>8}{'mode':>7}")
    for k in NAMES:
        m, sd, q, mode, _ = post(sel[k])
        A(f"{k:<18}{sel[k].simulation_id.nunique():>6}{m:>8.3f}{sd:>7.3f}" + "".join(f"{x:>8.2f}" for x in q) + f"{mode:>7.1f}")

P = {k: per_set(v[0]) for k, v in L.items()}
A("\n2. Run by run (each MC set's own KS-weighted mean age, and its lowest-KS age)")
A(f"{'pair':<30}{'sets':>6}{'median diff':>13}{'mean diff':>11}{'sd of diff':>12}{'|d|<0.25':>10}{'|d|<0.5':>9}{'corr':>7}")
for a, b in (("W v1.4 Z-matched", "DSEP3 physics on"), ("W v1.4 Z-matched", "OG"), ("DSEP3 physics on", "OG")):
    for kind in ("wmean", "best"):
        j = P[a][[kind]].join(P[b][[kind]], lsuffix="_a", rsuffix="_b", how="inner")
        d = j[f"{kind}_a"] - j[f"{kind}_b"]
        r = np.corrcoef(j[f"{kind}_a"], j[f"{kind}_b"])[0, 1]
        A(f"{f'{a} - {b} ({kind})':<30}{len(d):>6}{d.median():>+13.3f}{d.mean():>+11.3f}{d.std():>12.3f}"
          f"{(d.abs() < 0.25).mean():>10.0%}{(d.abs() < 0.5).mean():>9.0%}{r:>7.2f}")

j = P["W v1.4 Z-matched"][["wmean"]].join(P["DSEP3 physics on"][["wmean"]], lsuffix="_w", rsuffix="_d", how="inner")
diff = (j.wmean_w - j.wmean_d).rename("d")
X = iv.loc[diff.index]
A(f"\n3. What drives W v1.4 - DSEP3 (per-set weighted mean age), Spearman over {len(diff)} sets "
  f"(|rho| > {2 / np.sqrt(len(diff)):.2f} marked *):")
rows = []
for c in X.columns:
    if X[c].nunique() > 2:
        rho = spearmanr(X[c], diff).correlation
        rows.append((c, rho))
for c, rho in sorted(rows, key=lambda x: -abs(x[1]))[:10]:
    A(f"   {c:<32}{rho:+.2f}{' *' if abs(rho) > 2 / np.sqrt(len(diff)) else ''}")

A("\n4. Fit quality per set (lowest KS before the cut), sets fitted by both:")
for b in ("DSEP3 physics on", "OG"):
    k = pd.concat([L["W v1.4 Z-matched"][1].rename("w"), L[b][1].rename("o")], axis=1, join="inner")
    A(f"   W v1.4 vs {b}: median lowest KS {k.w.median():.4f} vs {k.o.median():.4f}; W v1.4 fits better in "
      f"{(k.w < k.o).mean():.0%} of {len(k)} sets")
A("\n5. By surface boundary condition (KTTAU), W v1.4 - DSEP3 per set and posterior on the sets all three keep:")
kt = iv.KTTAU
for v in sorted(kt.loc[diff.index].unique()):
    dd = diff[kt.loc[diff.index] == v]
    sel = [r for r in common if kt.loc[r] == v]
    ms = {k: post(L[k][0][L[k][0].simulation_id.isin(sel)])[:2] for k in NAMES} if sel else {}
    A(f"   KTTAU {int(v)}: {len(dd)} sets; per-set W v1.4 - DSEP3 median {dd.median():+.2f}, sd {dd.std():.2f} Gyr; more than 1 Gyr older "
      f"in W v1.4: {(dd > 1).sum()}" + ("; posterior (" + str(len(sel)) + " common sets): " + ", ".join(f"{k} {m:.2f} +/- {sd:.2f}" for k, (m, sd) in ms.items()) if sel else ""))
txt = "\n".join(out)
open(os.path.join(HERE, "age_breakdown.txt"), "w").write(txt + "\n")
print(txt)

# figure: posteriors on common sets, per-set W v1.4 vs DSEP3 ages, per-set differences
INK, MUTED, GRID, BG = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": BG, "axes.facecolor": BG})
fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
for k in NAMES:
    m, sd, q, mode, p = post(L[k][0][L[k][0].simulation_id.isin(common)])
    ax[0].step(ctr, p, where="mid", color=COL[k], lw=2, label=f"{k}: {m:.2f} ± {sd:.2f} Gyr")
ax[0].set_ylim(0, 0.115)
ax[0].set_xlim(5.5, 11.5)
ax[0].set_xlabel("age (Gyr)")
ax[0].set_ylabel("posterior probability per 0.1 Gyr")
ax[0].set_title(f"Posterior age, {len(common)} MC sets kept by all three", loc="left", fontsize=10)
ax[0].legend(frameon=False, fontsize=8.5, loc="upper left")
ax[1].plot([5, 12], [5, 12], color=MUTED, lw=1, ls="--")
KTC = {0: "#2a78d6", 5: "#1baf7a"}          # KTTAU 1 excluded (W v1.4 not trained on it)
for v, c in KTC.items():
    m_ = iv.KTTAU.loc[j.index] == v
    ax[1].plot(j.wmean_d[m_], j.wmean_w[m_], "o", ms=8, color=c, mec=BG, mew=1.5, alpha=0.9, label=f"KTTAU {v} ({m_.sum()} sets)")
ax[1].legend(frameon=False, fontsize=8.5, loc="upper left")
ax[1].set_xlabel("DSEP3 physics on: set's KS-weighted age (Gyr)")
ax[1].set_ylabel("W v1.4: same set (Gyr)")
ax[1].set_title(f"Same MC set, two models ({len(j)} sets; dashed = equal)", loc="left", fontsize=10)
ax[1].set_xlim(6, 12); ax[1].set_ylim(6, 12)
ktd = iv.KTTAU.loc[diff.index]
ax[2].hist([diff[ktd == v] for v in KTC], bins=np.arange(-3, 3.51, 0.25), color=list(KTC.values()), edgecolor=BG, lw=1.5,
           stacked=True, label=[f"KTTAU {v}" for v in KTC])
ax[2].legend(frameon=False, fontsize=8.5, loc="upper right")
ax[2].axvline(0, color=MUTED, lw=1, ls="--")
ax[2].axvline(diff.median(), color=INK, lw=1.5)
ax[2].annotate(f"median {diff.median():+.2f} Gyr\nsd {diff.std():.2f} Gyr", (diff.median(), 0), xytext=(-75, 120),
               textcoords="offset points", fontsize=8.5, color=INK)
ax[2].set_xlabel("W v1.4 − DSEP3 physics on, per MC set (Gyr)")
ax[2].set_ylabel("MC sets")
ax[2].set_title("Per-set age difference", loc="left", fontsize=10)
for a in ax:
    a.grid(axis="y", color=GRID, lw=0.8)
fig.suptitle(f"NGC 6791 sCMD ages: W v1.4, Z matched ({n_w14} runs) vs DSEP3 with physics on vs OG", x=0.01, ha="left",
             fontsize=11, color=INK)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(os.path.join(HERE, "age_breakdown.png"), dpi=150)
print(os.path.join(HERE, "age_breakdown.png"))
