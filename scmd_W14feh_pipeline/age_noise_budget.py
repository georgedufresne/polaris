"""NGC 6791 sCMD age: where the posterior width comes from, and how much of the W v1.4 - DSEP3 per-set scatter is
isochrone (emulator) error rather than sCMD Monte Carlo noise.

Same cut (0.0468), weights (KDE of the original bootstrap null), one best DM/E per (MC set, age), KTTAU 0/5.

1. Law of total variance per leg: posterior variance = weighted mean of within-set variances (the fit: age traded
   against DM/E(B-V) within the KS acceptance) + weighted variance of the set means (between sets). The between-set part
   is split by a weighted linear fit of set mean age on the var-file inputs (physics) and its residual.
2. Repeatability of the sCMD pipeline: the reruns of 2026-11-22 recomputed whole runs with the same isochrones (new
   synthetic stars). The first attempts' partial tables (2d_ks_noco_partial, 2d_ks_W14_partial) are compared with the
   reruns on the ages present in both: per-set KS-weighted mean age, attempt 1 - attempt 2. sd / sqrt(2) = sCMD noise
   of one set's age.
3. Emulator error: per-set W14feh - New sd (equal weight) = sqrt(2 s_sCMD^2 + s_iso^2) if the sCMD noise is
   independent between legs, giving s_iso, the isochrone (emulator) part.
Writes age_noise_budget.txt here.
"""
import glob
import os

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

S = "/dartfs-hpc/scratch/f007058"
HERE = os.path.dirname(os.path.abspath(__file__))
CUT = 0.0468
RNG = np.random.default_rng(1123)
kde = gaussian_kde(pd.read_csv(f"{HERE}/inputs/bootstrap_2d_ks_results.csv").ks_statistic.values, bw_method="scott")
IV = pd.read_csv(f"{S}/analysis_W/initial_vars.csv").rename(columns={"file_id": "simulation_id"}).set_index("simulation_id")
KEEP = set(IV.index[IV.KTTAU.isin([0, 5])])
COLS = [n for n in IV.columns if IV[n].std() > 0]


def rows(d):
    d = d[d.simulation_id.isin(KEEP) & (d.ks_statistic <= CUT)].copy()
    d["w"] = kde(d.ks_statistic.values)
    return d.loc[d.groupby(["simulation_id", "age_gyr"]).ks_statistic.idxmin()].reset_index(drop=True)


def set_stats(d):
    g = d.groupby("simulation_id")
    m = pd.Series({k: np.average(x.age_gyr, weights=x.w) for k, x in g})
    v = pd.Series({k: np.average((x.age_gyr - m[k]) ** 2, weights=x.w) for k, x in g})
    return pd.DataFrame({"m": m, "v": v, "W": g.w.sum()})


LEGS = {"OG": f"{HERE}/inputs/filtered_ks_scmd_OG.csv", "New": f"{S}/analysis_noco/filtered_ks_scmd.csv",
        "W14": f"{S}/analysis_W14/filtered_ks_scmd.csv", "W14feh": f"{S}/analysis_W14feh/filtered_ks_scmd.csv"}
ST = {k: set_stats(rows(pd.read_csv(p))) for k, p in LEGS.items()}
out = []
A = out.append
A("NGC 6791 sCMD age noise budget (KTTAU 0/5; weights and cut as scmd_6way_compare.py)\n")
A("1. Posterior age variance = within-set (fit) + between-set (inputs explained + residual), each as an sd in Gyr")
A(f"  {'leg':<8}{'sets':>6}{'total sd':>10}{'within':>9}{'between':>9}{'  inputs':>9}{'residual':>10}{'  R^2 inputs':>13}")
for k, t in ST.items():
    p = t.W / t.W.sum()
    mu = (p * t.m).sum()
    within = (p * t.v).sum()
    between = (p * (t.m - mu) ** 2).sum()
    X = np.column_stack([np.ones(len(t))] + [IV.loc[t.index, n].values for n in COLS])
    sw = np.sqrt(p.values)
    beta, *_ = np.linalg.lstsq(X * sw[:, None], t.m.values * sw, rcond=None)
    yh = X @ beta
    expl = (p * (yh - mu) ** 2).sum()
    res = (p * (t.m - yh) ** 2).sum()
    A(f"  {k:<8}{len(t):>6}{np.sqrt(within + between):>10.3f}{np.sqrt(within):>9.3f}{np.sqrt(between):>9.3f}"
      f"{np.sqrt(expl):>9.3f}{np.sqrt(res):>10.3f}{expl / between:>13.2f}")
A("  (variances add: total^2 = within^2 + between^2; between^2 = inputs^2 + residual^2)\n")

A("2. sCMD repeatability: same isochrones, first attempt (partial table) vs rerun, on the ages present in both")
rep = []
for part, new, leg in [("2d_ks_noco_partial", "2d_ks_noco", "New"), ("2d_ks_W14_partial", "2d_ks_W14", "W14")]:
    for f in sorted(glob.glob(f"{S}/{part}/2d_ks_results_*.csv")):
        g = f"{S}/{new}/{os.path.basename(f)}"
        if not os.path.exists(g):
            continue
        a, b = pd.read_csv(f), pd.read_csv(g)
        ages = set(a.age_gyr.round(1)) & set(b.age_gyr.round(1))
        if len(ages) < 60:
            continue
        a, b = a[a.age_gyr.round(1).isin(ages)], b[b.age_gyr.round(1).isin(ages)]
        ra, rb = rows(a), rows(b)
        if ra.empty or rb.empty:
            continue
        sa, sb = set_stats(ra), set_stats(rb)
        r = int(a.simulation_id.iloc[0])
        mina = a.groupby(a.age_gyr.round(1)).ks_statistic.min()
        minb = b.groupby(b.age_gyr.round(1)).ks_statistic.min()
        rep.append((leg, r, len(ages), sa.m.iloc[0], sb.m.iloc[0], (mina - minb).std(), a.ks_statistic.min(), b.ks_statistic.min()))
rep = pd.DataFrame(rep, columns=["leg", "run", "ages", "age1", "age2", "sd_minKS_diff", "best1", "best2"])
rep["d"] = rep.age1 - rep.age2
for _, x in rep.iterrows():
    A(f"  {x.leg:<4} run {x.run}: {x.ages} ages, set age {x.age1:.3f} vs {x.age2:.3f} (diff {x.d:+.3f}); best KS {x.best1:.4f} vs "
      f"{x.best2:.4f}; per-age min-KS diff sd {x.sd_minKS_diff:.4f}")
s_rep = rep.d.std(ddof=1) / np.sqrt(2)
bs = np.array([rep.d.iloc[RNG.integers(0, len(rep), len(rep))].std(ddof=1) for _ in range(5000)]) / np.sqrt(2)
mad = 1.4826 * np.median(np.abs(rep.d - rep.d.median())) / np.sqrt(2)
A(f"  {len(rep)} runs kept by the cut in both attempts: sCMD noise of one set's age = {s_rep:.3f} Gyr "
  f"(bootstrap 68% {np.percentile(bs, 16):.3f}-{np.percentile(bs, 84):.3f}; robust MAD estimate {mad:.3f})\n")

A("3. W14feh - New per set (equal weight)")
j = ST["W14feh"].join(ST["New"], lsuffix="_w", rsuffix="_d", how="inner")
d = j.m_w - j.m_d
sd_d = d.std(ddof=1)
mad_d = 1.4826 * np.median(np.abs(d - d.median()))
A(f"  {len(j)} sets: mean {d.mean():+.3f}, median {d.median():+.3f}, sd {sd_d:.3f}, robust sd (MAD) {mad_d:.3f}")
for lab, sdd, sr in [("sd", sd_d, s_rep), ("robust", mad_d, mad)]:
    iso2 = sdd ** 2 - 2 * sr ** 2
    A(f"  {lab}: sCMD noise in the difference sqrt(2) x {sr:.3f} = {np.sqrt(2) * sr:.3f}; "
      f"isochrone (emulator) part = {np.sqrt(iso2) if iso2 > 0 else 0:.3f} Gyr per set")
A("  W14 - New for comparison: " + (lambda dd: f"mean {dd.mean():+.3f}, median {dd.median():+.3f}, sd {dd.std(ddof=1):.3f}")(
    (lambda jj: jj.m_w - jj.m_d)(ST["W14"].join(ST["New"], lsuffix="_w", rsuffix="_d", how="inner"))))
open(f"{HERE}/age_noise_budget.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
