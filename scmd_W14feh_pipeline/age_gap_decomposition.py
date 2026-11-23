"""Pooled posterior age difference (a - b) on the MC sets both legs keep, split into a within-set shift (sum of
p_a,s * (m_a,s - m_b,s)) and a set-reweighting term (sum of (p_a,s - p_b,s) * m_b,s). m = each set's KS-weighted
mean age, p = its share of the leg's total weight. Same cut, weights and KTTAU 0/5 restriction as scmd_6way_compare.py.
Writes age_gap_decomposition.txt here."""
import glob, numpy as np, pandas as pd
from scipy.stats import gaussian_kde
S = "/dartfs-hpc/scratch/f007058"; CUT = 0.0468
kde = gaussian_kde(pd.read_csv(f"{S}/scmd_W14feh_pipeline/inputs/bootstrap_2d_ks_results.csv").ks_statistic.values, bw_method="scott")
iv = pd.read_csv(f"{S}/analysis_W/initial_vars.csv").rename(columns={"file_id": "simulation_id"}).set_index("simulation_id")
KEEP = set(iv.index[iv.KTTAU.isin([0, 5])])
def prep(d):
    d = d[d.simulation_id.isin(KEEP) & (d.ks_statistic <= CUT)].copy(); d["w"] = kde(d.ks_statistic.values)
    d = d.loc[d.groupby(["simulation_id", "age_gyr"]).ks_statistic.idxmin()]
    g = d.groupby("simulation_id")
    return pd.DataFrame({"m": pd.Series({k: np.average(x.age_gyr, weights=x.w) for k, x in g}), "W": g.w.sum(), "n": g.size()})
raw = lambda k: pd.concat([pd.read_csv(f) for f in glob.glob(f"{S}/{k}/2d_ks_results_*.csv")], ignore_index=True)
L = {"New": prep(pd.read_csv(f"{S}/analysis_noco/filtered_ks_scmd.csv")), "W14feh": prep(raw("2d_ks_W14feh")),
     "W14": prep(pd.read_csv(f"{S}/analysis_W14/filtered_ks_scmd.csv")), "OG": prep(pd.read_csv(f"{S}/scmd_W14feh_pipeline/inputs/filtered_ks_scmd_OG.csv"))}
out = []
for a, b in [("W14feh", "New"), ("W14", "New"), ("W14feh", "OG"), ("New", "OG")]:
    j = L[a].join(L[b], lsuffix="a", rsuffix="b", how="inner")
    pa, pb = j.Wa / j.Wa.sum(), j.Wb / j.Wb.sum()
    pooled = (pa * j.ma).sum() - (pb * j.mb).sum()
    within = (pa * (j.ma - j.mb)).sum(); rew = ((pa - pb) * j.mb).sum()
    eq = (j.ma - j.mb).mean()
    c = np.corrcoef(np.log(j.Wa / j.Wb), j.mb)[0, 1]
    out.append(f"{a}-{b}: n={len(j)} pooled {pooled:+.3f} = within-set {within:+.3f} + reweighting {rew:+.3f}; equal-weight mean {eq:+.3f}; "
          f"rows/set {j.na.mean():.1f} vs {j.nb.mean():.1f}; corr(log Wa/Wb, age_b) {c:+.2f}")
open(f"{S}/scmd_W14feh_pipeline/age_gap_decomposition.txt", "w").write("\n".join(out) + "\n")
print("\n".join(out))
