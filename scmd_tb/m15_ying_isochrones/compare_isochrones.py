#!/usr/bin/env python3
"""Ying et al. (2025)'s published M15 Monte Carlo isochrones (Zenodo 10.5281/zenodo.13998637, M15_MC_iso.h5)
against our DSEP4 rebuild (iso_dsep4, MIST ACS/WFC BCs at [Fe/H] -2.25) for the same runs.

1. Inputs: his 'var' table vs our var-derived inputs (posterior_weights.csv columns) for each run.
2. Same run, same age: turnoff (bluest point) F606W and colour, the SGB point 0.05 mag redder than the turnoff
   (Ying's V_SGB definition), the main-sequence colour at M_F606W 5.8 and 4.5, the RGB colour 1.5 mag above the SGB
   point. Differences are his minus ours.
3. Equivalent age: his isochrone at 13.4 Gyr, within +/-2 mag of its SGB point (the fit window), matched by each of
   our ages with a free F606W shift and a free colour shift (as distance and reddening are free in the fit),
   minimising the mean squared distance in (6 x colour, F606W). The best-matching age of ours is reported.

His isochrones are read with h5py (pandas 'fixed' format; data only). Output: logs/compare_isochrones.out, compare_isochrones.csv."""
import json
import sys

import h5py
import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.spatial import cKDTree

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
H = f"{B}/m15_ying_isochrones"
YH5 = "/dartfs-hpc/scratch/f007058/M15/zenodo_ying2025/M15_MC_iso.h5"
sys.path.insert(0, f"{B}/m15_final/code")
from isochronetoolbox.scmd.model import BCTable, read_isochrones  # noqa: E402

K = 6.0                                    # colour weight in the curve match (about the main-sequence slope)
lim = {int(k): v for k, v in json.load(open("/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json")).items()}
bc = BCTable(f"{B}/m15_grid_d4/inputs/bc_mist_acs_wfc_fehm225.npz")
w = pd.read_csv(f"{B}/m15_final_refined/posterior_scaled/posterior_weights.csv")
ours_var = w.drop_duplicates("simulation_id").set_index("simulation_id")
wt = w.groupby("simulation_id").P_KS_normalized.sum().sort_values(ascending=False)
RUNS = sorted(set([10357, 10392, 10245] + list(wt.index[:27])))       # key runs + the 27 most-weighted in production


YF = h5py.File(YH5, "r")


def read_df(path):
    """A pandas 'fixed'-format frame (block0_items / block0_values) read with h5py, without PyTables."""
    g = YF[path]
    return pd.DataFrame(g["block0_values"][()], columns=[x.decode() for x in g["block0_items"][()]])


def ying(run, age):
    d = read_df(f"{run}/isochrones/cmd/{int(round(age * 1000))}").sort_values("EEP")
    return d.F606W.values, d.F814W.values                           # 'F814W' column holds F606W-F814W


def ours_all(run):
    iso = read_isochrones(f"/dartfs-hpc/scratch/f007058/M15/iso_dsep4/{run}/isos/isochrones_example.txt")
    m = bc.absolute(iso)
    iso = iso.assign(v=m["F606W"], c=m["F606W"] - m["F814W"])
    return {round(a / 1e9, 1): (g.sort_values("EEP").v.values, g.sort_values("EEP").c.values)
            for a, g in iso.groupby("isochrone_age_yr")}


def features(v, c):
    up = np.where(v < 5.5)[0]
    i_to = up[np.argmin(c[up])]
    post = np.where((np.arange(len(c)) > i_to) & (c >= c[i_to] + 0.05))[0]
    i_sgb = post[0]
    ms = np.where(np.arange(len(c)) < i_to)[0]
    rgb = np.where(np.arange(len(c)) > i_sgb)[0]
    vs = v[i_sgb]
    return dict(v_to=v[i_to], c_to=c[i_to], v_sgb=vs,
                c_ms58=np.interp(5.8, v[ms][::-1], c[ms][::-1]), c_ms45=np.interp(4.5, v[ms][::-1], c[ms][::-1]),
                c_rgb=np.interp(vs - 1.5, v[rgb][::-1], c[rgb][::-1]) if len(rgb) and v[rgb].min() < vs - 1.5 else np.nan)


def dense(v, c, step=0.002):
    s = np.concatenate([[0], np.cumsum(np.hypot(K * np.diff(c), np.diff(v)))])
    t = np.arange(0, s[-1], step)
    return np.column_stack([K * np.interp(t, s, c), np.interp(t, s, v)])


def match(target, ours_curve, x0):
    tree = cKDTree(dense(*ours_curve))
    f = lambda p: np.mean(tree.query(target + [K * p[1], p[0]])[0] ** 2)
    best = minimize(f, x0, method="Nelder-Mead", options=dict(xatol=1e-4, fatol=1e-10))
    return best.fun, best.x


rows = []
VARS = [("FeH", "FeH"), ("Yprim", "Y"), ("alphafe", "alphafe"), ("CMIXLA", "CMIXLA"), ("FGRY", "FGRY"), ("FGRZ", "FGRZ"),
        ("ALPHAE", "ALPHAE"), ("ALPHAC", "ALPHAC")]
for run in RUNS:
    yv = read_df(f"{run}/var").iloc[0]
    dv = {a: yv[a] - ours_var.loc[run, b] for a, b in VARS}
    ours = ours_all(run)
    r = dict(simulation_id=run, FeH=yv['FeH'], CMIXLA=yv['CMIXLA'], ALPHAC_ying=yv['ALPHAC'], alphafe=yv['alphafe'],
             FGRY=yv['FGRY'], Yprim=yv['Yprim'], **{f"d_{a}": x for a, x in dv.items()})
    for age in (13.4, 14.2):
        fy, fo = features(*ying(run, age)), features(*ours[age])
        r.update({f"{k}_{age}": fy[k] - fo[k] for k in fy})
        r[f"ying_v_sgb_minus_v_to_{age}"] = fy["v_sgb"] - fy["v_to"]
        r[f"ours_v_sgb_minus_v_to_{age}"] = fo["v_sgb"] - fo["v_to"]
    v, c = ying(run, 13.4)
    vs = features(v, c)["v_sgb"]
    k = (v > vs - 2) & (v < vs + 2)
    target = np.column_stack([K * c[k], v[k]])
    fy = features(v, c)
    ok = [a for a in sorted(ours) if 11.0 <= a <= 15.8 and a < lim[run] - 1e-6]

    def fit_age(a):
        fo = features(*ours[a])          # start from the shift that lines up the SGB point and the turnoff colour
        return (a, *match(target, ours[a], [fy["v_sgb"] - fo["v_sgb"], fy["c_to"] - fo["c_to"]]))
    coarse = [fit_age(a) for a in ok[::2]]
    a0 = min(coarse, key=lambda x: x[1])[0]
    fits = coarse + [fit_age(a) for a in ok if abs(a - a0) < 0.25 and a not in ok[::2]]
    a, cost, (dm, dc) = min(fits, key=lambda x: x[1])
    r.update(equiv_age_for_ying_13_4=a, equiv_rms=np.sqrt(cost), equiv_dv=dm, equiv_dc=dc)
    rows.append(r)
    print(f"run {run}: inputs max |diff| {max(abs(x) for x in dv.values()):.2g} | at 13.4 Gyr his-ours: TO F606W "
          f"{r['v_to_13.4']:+.3f} colour {r['c_to_13.4']:+.4f}, SGB F606W {r['v_sgb_13.4']:+.3f}, MS colour at 4.5 "
          f"{r['c_ms45_13.4']:+.4f}, RGB colour {r['c_rgb_13.4']:+.4f} | our age matching his 13.4: {a:.1f} "
          f"(shift dF606W {dm:+.3f}, dcolour {dc:+.4f}, rms {np.sqrt(cost):.4f})", flush=True)
df = pd.DataFrame(rows)
df.to_csv(f"{H}/compare_isochrones.csv", index=False)
print(f"\n{len(df)} runs. Median (16-84%) of his minus ours, same run and age:")
for col in ["v_to", "c_to", "v_sgb", "c_ms58", "c_ms45", "c_rgb"]:
    for age in (13.4, 14.2):
        x = df[f"{col}_{age}"]
        print(f"  {col:7s} {age}: {x.median():+.4f} ({x.quantile(.16):+.4f} to {x.quantile(.84):+.4f})")
for age in (13.4, 14.2):
    print(f"  SGB-TO F606W gap at {age}: his {df[f'ying_v_sgb_minus_v_to_{age}'].median():.3f}, ours {df[f'ours_v_sgb_minus_v_to_{age}'].median():.3f}")
e = df.equiv_age_for_ying_13_4
print(f"our age that best matches his 13.4 Gyr isochrone (free shifts): median {e.median():.2f}, 16-84% "
      f"{e.quantile(.16):.2f}-{e.quantile(.84):.2f}; shifts median dF606W {df.equiv_dv.median():+.3f}, "
      f"dcolour {df.equiv_dc.median():+.4f}; rms {df.equiv_rms.median():.4f}")
print("input differences (his - ours), max |diff| over runs: "
      + ", ".join(f"{a} {df[f'd_{a}'].abs().max():.3g}" for a, _ in VARS))

print("Spearman correlation with the inputs, of the turnoff colour difference (13.4 Gyr) and of the equivalent age: "
      + "; ".join(f"{c} {df['c_to_13.4'].corr(df[c], method='spearman'):+.2f} / {e.corr(df[c], method='spearman'):+.2f}"
                  for c in ["FeH", "CMIXLA", "ALPHAC_ying", "alphafe", "FGRY", "Yprim"]))
