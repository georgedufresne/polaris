"""NGC 6791 sCMD: every parameter, paired leg against leg, and what drives the remaining W v1.4 - DSEP3 age gap.

Same cut (0.0468), weights (KDE of the original bootstrap null at each KS value), one best DM/E per (MC set, age), and
KTTAU 0/5 restriction as scmd_6way_compare.py.

1. Paired pooled difference on the MC sets both legs keep, for age, DM, E(B-V) (fitted per row) and [Fe/H], Y,
   [alpha/Fe], mixing length (var-file inputs, fixed per set). Bootstrap SE over sets. Each is split into
   within-set = sum p_a,s (m_a,s - m_b,s) and reweighting = sum (p_a,s - p_b,s) m_b,s, where m is the set's KS-weighted
   mean and p its share of the leg's weight. The inputs have no within-set part, so their difference is pure reweighting.
2. Within-set (W14feh - New): per-set age shift against the per-set DM and E(B-V) shifts (does W trade age for
   distance/reddening?), and a linear fit dage = a + b dDM + c dE.
3. Reweighting (W14feh vs New): Spearman of log(W_W14feh / W_New) per set with each input and with DSEP3's per-set age;
   the reweighting term split into the part predicted by a linear fit of DSEP3's per-set age on all var-file inputs
   and the residual. Bootstrap SE of the within-set and reweighting age terms for each pair.
0. Age posterior per leg (own sets; sets kept by all four legs): weighted mean +/- sd, median, 16-84 % and 2.5-97.5 %
   ranges, from the same 0.1 Gyr histogram as age_breakdown_W14feh.py. Per pair: pooled mean and median differences, and
   the per-set (equal weight) mean and median differences, each with bootstrap SE over sets.
Writes param_gap_breakdown.txt here.
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, spearmanr

S = "/dartfs-hpc/scratch/f007058"
HERE = os.path.dirname(os.path.abspath(__file__))
CUT, NB = 0.0468, 2000
RNG = np.random.default_rng(1123)
kde = gaussian_kde(pd.read_csv(f"{HERE}/inputs/bootstrap_2d_ks_results.csv").ks_statistic.values, bw_method="scott")
IV = pd.read_csv(f"{S}/analysis_W/initial_vars.csv").rename(columns={"file_id": "simulation_id"}).set_index("simulation_id")
KEEP = set(IV.index[IV.KTTAU.isin([0, 5])])
FIT = [("age (Gyr)", "age_gyr"), ("DM", "distance_modulus"), ("E(B-V)", "extinction_ebv")]
INP = [("[Fe/H]", "FeH"), ("Y", "He Abundance"), ("[alpha/Fe]", "alphafe"), ("mixing length", "CMIXLA")]


ROWS = {}
BINS = np.arange(3.95, 12.15, 0.1)
CTR = (BINS[:-1] + BINS[1:]) / 2


def post(age, w):
    h, _ = np.histogram(age, bins=BINS, weights=w)
    p = h / h.sum()
    m = (CTR * p).sum()
    c = np.cumsum(p)
    q = [np.interp(x, c, CTR) for x in (0.5, 0.16, 0.84, 0.025, 0.975)]
    return m, np.sqrt(((CTR - m) ** 2 * p).sum()), q


def per_set(path, name):
    d = pd.read_csv(path)
    d = d[d.simulation_id.isin(KEEP) & (d.ks_statistic <= CUT)].copy()
    d["w"] = kde(d.ks_statistic.values)
    d = d.loc[d.groupby(["simulation_id", "age_gyr"]).ks_statistic.idxmin()]
    ROWS[name] = d[["simulation_id", "age_gyr", "w"]].reset_index(drop=True)
    g = d.groupby("simulation_id")
    out = {c: pd.Series({k: np.average(x[c], weights=x.w) for k, x in g}) for _, c in FIT}
    out["W"] = g.w.sum()
    out = pd.DataFrame(out)
    for _, c in INP:
        out[c] = IV.loc[out.index, c].values
    return out


L = {"OG": per_set(f"{HERE}/inputs/filtered_ks_scmd_OG.csv", "OG"), "New": per_set(f"{S}/analysis_noco/filtered_ks_scmd.csv", "New"),
     "W14": per_set(f"{S}/analysis_W14/filtered_ks_scmd.csv", "W14"),
     "W14feh": per_set(f"{S}/analysis_W14feh/filtered_ks_scmd.csv", "W14feh")}
out = []
A = out.append


def split(a, b, c):
    pa, pb = a.W / a.W.sum(), b.W / b.W.sum()
    return (pa * a[c]).sum() - (pb * b[c]).sum(), (pa * (a[c] - b[c])).sum(), ((pa - pb) * b[c]).sum()


A("NGC 6791 sCMD, KTTAU 0/5: every parameter paired on common MC sets (pooled = within-set + reweighting; bootstrap SE)")
A("Legs: OG (physics off), New (DSEP3 physics on), W14 (W v1.4, var [Fe/H]), W14feh (W v1.4, [Fe/H] for our Z)\n")
all4 = set.intersection(*[set(v.index) for v in L.values()])
for title, keep in [("each leg's own sets", None), (f"sets kept by all four legs ({len(all4)})", all4)]:
    A(f"Age posterior (Gyr), {title}:")
    A(f"  {'leg':<8}{'sets':>6}{'mean':>8}{'sd':>7}{'median':>8}{'16-84 %':>15}{'2.5-97.5 %':>15}{'median-mean':>13}")
    for k, r in ROWS.items():
        r = r if keep is None else r[r.simulation_id.isin(keep)]
        m, sd, q = post(r.age_gyr, r.w)
        A(f"  {k:<8}{r.simulation_id.nunique():>6}{m:>8.3f}{sd:>7.3f}{q[0]:>8.3f}{f'{q[1]:.2f}-{q[2]:.2f}':>15}{f'{q[3]:.2f}-{q[4]:.2f}':>15}{q[0] - m:>+13.3f}")
    A("")
for an, bn in [("W14feh", "New"), ("W14", "New"), ("W14feh", "W14"), ("New", "OG"), ("W14feh", "OG")]:
    common = L[an].index.intersection(L[bn].index)
    a, b = L[an].loc[common], L[bn].loc[common]
    A(f"{an} - {bn}  ({len(common)} common sets)")
    A(f"  {'parameter':<15}{'pooled':>10}{'SE':>9}{'z':>7}{'within-set':>12}{'reweighting':>13}")
    for lab, c in FIT + INP:
        pooled, within, rew = split(a, b, c)
        bs = np.empty(NB)
        for i in range(NB):
            k = RNG.integers(0, len(common), len(common))
            bs[i] = split(a.iloc[k], b.iloc[k], c)[0]
        se = bs.std(ddof=1)
        A(f"  {lab:<15}{pooled:>+10.4f}{se:>9.4f}{pooled / se:>7.2f}{within:>+12.4f}{rew:>+13.4f}")
    bw = np.array([split(a.iloc[k], b.iloc[k], "age_gyr")[1:] for k in
                   (RNG.integers(0, len(common), len(common)) for _ in range(NB))])
    A(f"  age terms, bootstrap SE: within-set +/- {bw[:, 0].std(ddof=1):.4f}, reweighting +/- {bw[:, 1].std(ddof=1):.4f}")
    ra, rb = ROWS[an][ROWS[an].simulation_id.isin(common)], ROWS[bn][ROWS[bn].simulation_id.isin(common)]
    ga, gb = ra.groupby("simulation_id").indices, rb.groupby("simulation_id").indices
    sets = np.asarray(common)
    d_set = (a.age_gyr - b.age_gyr).values

    def stats(pick):
        ia = np.concatenate([ga[sets[j]] for j in pick]); ib = np.concatenate([gb[sets[j]] for j in pick])
        pa_, pb_ = post(ra.age_gyr.values[ia], ra.w.values[ia]), post(rb.age_gyr.values[ib], rb.w.values[ib])
        return pa_[0] - pb_[0], pa_[2][0] - pb_[2][0], d_set[pick].mean(), np.median(d_set[pick])
    obs = stats(np.arange(len(sets)))
    bs = np.array([stats(RNG.integers(0, len(sets), len(sets))) for _ in range(NB)])
    for lab, o, x in zip(["posterior mean diff", "posterior median diff", "per-set mean diff", "per-set median diff"], obs, bs.T):
        lo, hi = np.percentile(x, [2.5, 97.5])
        A(f"  age {lab:<22}{o:>+8.3f} +/- {x.std(ddof=1):.3f}   95% CI [{lo:+.3f}, {hi:+.3f}]")
    A("")

common = L["W14feh"].index.intersection(L["New"].index)
a, b = L["W14feh"].loc[common], L["New"].loc[common]
dd = pd.DataFrame({c: a[c] - b[c] for _, c in FIT})
A("Within-set, W14feh - New (per-set KS-weighted means, equal weight per set):")
for lab, c in FIT:
    A(f"  {lab:<10} mean {dd[c].mean():+.4f}  median {dd[c].median():+.4f}  sd {dd[c].std():.4f}")
A(f"  corr(d age, d DM) {dd.age_gyr.corr(dd.distance_modulus):+.2f}   corr(d age, d E(B-V)) {dd.age_gyr.corr(dd.extinction_ebv):+.2f}   "
  f"corr(d DM, d E(B-V)) {dd.distance_modulus.corr(dd.extinction_ebv):+.2f}")
X = np.column_stack([np.ones(len(dd)), dd.distance_modulus, dd.extinction_ebv])
coef, *_ = np.linalg.lstsq(X, dd.age_gyr.values, rcond=None)
A(f"  fit: d age = {coef[0]:+.4f} {coef[1]:+.3f} x d DM {coef[2]:+.3f} x d E(B-V)   (intercept = age shift at equal DM and E(B-V))")
A(f"  share of the mean age shift carried by DM/E shifts: {(coef[1] * dd.distance_modulus.mean() + coef[2] * dd.extinction_ebv.mean()):+.4f} Gyr\n")

lr = np.log(a.W / b.W)
A("Reweighting, W14feh vs New: Spearman of log(weight ratio) per set with")
rows = [(lab, spearmanr(lr, b[c]).statistic) for lab, c in INP] + [("DSEP3 per-set age", spearmanr(lr, b.age_gyr).statistic)]
rows += [(n, spearmanr(lr, IV.loc[common, n]).statistic) for n in IV.columns if n not in [c for _, c in INP] + ["KTTAU"]]
for n, r in sorted(rows, key=lambda t: -abs(t[1]))[:10]:
    A(f"  {n:<32}{r:+.2f}")
A("  and DSEP3 per-set age with the inputs: " + ", ".join(f"{lab} {spearmanr(b.age_gyr, b[c]).statistic:+.2f}" for lab, c in INP))

cols = [n for n in IV.columns if IV.loc[common, n].std() > 0]
Xi = np.column_stack([np.ones(len(common))] + [IV.loc[common, n].values for n in cols])
beta, *_ = np.linalg.lstsq(Xi, b.age_gyr.values, rcond=None)
yhat = Xi @ beta
pa, pb = a.W / a.W.sum(), b.W / b.W.sum()
r2 = 1 - ((b.age_gyr - yhat) ** 2).sum() / ((b.age_gyr - b.age_gyr.mean()) ** 2).sum()
A(f"  linear fit of DSEP3 per-set age on {len(cols)} inputs: R^2 {r2:.2f}")
A(f"  reweighting {((pa - pb) * b.age_gyr).sum():+.4f} = through the inputs {((pa - pb) * yhat).sum():+.4f} + residual {((pa - pb) * (b.age_gyr - yhat)).sum():+.4f}")
contrib = sorted(((n, beta[i + 1] * ((pa - pb) * IV.loc[common, n].values).sum(), beta[i + 1] * IV.loc[common, n].std(),
                  ((pa - pb) * IV.loc[common, n].values).sum() / IV.loc[common, n].std()) for i, n in enumerate(cols)), key=lambda t: -abs(t[1]))
A("  by input (age contribution = fit slope x weight shift; slope x 1 sd; weight shift of the input in sd):")
for n, c, sl, sh in contrib[:8]:
    A(f"    {n:<32}{c:+.4f} Gyr   (age per sd {sl:+.3f}, shift {sh:+.3f} sd)")
mlt = np.polyfit(b.CMIXLA, b.age_gyr, 1)
A(f"  mixing length alone: DSEP3 age slope {mlt[0]:+.2f} Gyr per unit; reweighting through it {((pa - pb) * np.polyval(mlt, b.CMIXLA)).sum():+.4f}")

open(f"{HERE}/param_gap_breakdown.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
