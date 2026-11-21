#!/usr/bin/env python3
"""Calibration-star weights for the M15 Monte Carlo runs (Ying et al. 2023, M92 paper, sections 2.1 and 3.1), and
the M15 age posteriors reweighted by them.

Stars (Ying et al. 2023, section 2.1): HIP 46120 M_F606W 5.7867 +/- 0.0026, F606W-F814W 0.566 +/- 0.002;
HIP 106924 M_F606W 6.0406 +/- 0.0037, F606W-F814W 0.601 +/- 0.005. Zero reddening, Gaia EDR3 parallaxes, "virtually
the identical composition as M92" ([Fe/H] about -2.3); the 2025 paper lists the same two stars for M15.

Recipe (Ying et al. 2023): for each age of each run, chi2 = sum over the two stars of the squared, uncertainty-
normalised F606W and colour differences between the isochrone and the star. Here that is the squared normalised
distance from the star to the isochrone curve (piecewise linear between EEPs; the colour difference at the star's
magnitude is printed as well). The run's weight is p = 1 - CDF_chi2(n=2)(chi2_min) = exp(-chi2_min / 2), with
chi2_min the minimum over the run's ages (8-16 Gyr, capped at the run's age limit). Our isochrones and bolometric
corrections exactly as in production (iso_dsep4 isochrones_example.txt, MIST ACS/WFC BCs at [Fe/H] -2.25).

The posterior weights (P_KS_normalized in each posterior_weights.csv) are multiplied by p(run) and renormalised.
Output: calib_chi2.csv (per run), logs/calib_weights.out."""
import json
import sys

import numpy as np
import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
H = f"{B}/m15_calib_stars"
sys.path.insert(0, f"{B}/m15_final/code")
from isochronetoolbox.scmd.model import BCTable, read_isochrones  # noqa: E402

STARS = {"HIP 46120": (5.7867, 0.0026, 0.566, 0.002), "HIP 106924": (6.0406, 0.0037, 0.601, 0.005)}
lim = {int(k): v for k, v in json.load(open("/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json")).items()}
runs = [int(x) for x in open(f"{B}/m15_final_refined/runs.txt").read().split()]
bc = BCTable(f"{B}/m15_grid_d4/inputs/bc_mist_acs_wfc_fehm225.npz")


def dist2(x, y, x0, y0):
    """Squared distance from (x0, y0) to the polyline (x, y)."""
    ax, ay, bx, by = x[:-1], y[:-1], x[1:], y[1:]
    dx, dy = bx - ax, by - ay
    t = np.clip(((x0 - ax) * dx + (y0 - ay) * dy) / np.maximum(dx * dx + dy * dy, 1e-30), 0, 1)
    return np.min((ax + t * dx - x0) ** 2 + (ay + t * dy - y0) ** 2)


rows = []
for run in runs:
    iso = read_isochrones(f"/dartfs-hpc/scratch/f007058/M15/iso_dsep4/{run}/isos/isochrones_example.txt")
    m = bc.absolute(iso)
    iso = iso.assign(v=m["F606W"], c=m["F606W"] - m["F814W"])
    best = None
    for age_yr, g in iso.groupby("isochrone_age_yr"):
        age = age_yr / 1e9
        if age < 8.0 - 1e-6 or age >= min(16.0 + 1e-6, lim[run] - 1e-6):
            continue
        g = g.sort_values("EEP")
        chi2, dc = 0.0, {}
        for name, (M, sM, C, sC) in STARS.items():
            chi2 += dist2(g.c.values / sC, g.v.values / sM, C / sC, M / sM)
            ms = g[(g.v > 4.5) & (g.v < 7.5)].sort_values("v")
            dc[name] = np.interp(M, ms.v.values, ms.c.values) - C      # isochrone colour at the star's F606W - star
        if best is None or chi2 < best[1]:
            best = (age, chi2, dc)
    age, chi2, dc = best
    rows.append(dict(simulation_id=run, chi2_min=chi2, age_at_min=age, p=np.exp(-chi2 / 2),
                     dc_46120=dc["HIP 46120"], dc_106924=dc["HIP 106924"]))
cal = pd.DataFrame(rows)
cal.to_csv(f"{H}/calib_chi2.csv", index=False)

print(f"{len(cal)} runs. chi2_min: median {cal.chi2_min.median():.1f}, 16-84% {cal.chi2_min.quantile(.16):.1f}-"
      f"{cal.chi2_min.quantile(.84):.1f}; chi2 < 2: {(cal.chi2_min < 2).mean():.0%} (M92 paper: about 20%); "
      f"p > 0.05: {(cal.p > 0.05).mean():.0%} (M92 paper: about 35%)")
print(f"isochrone colour minus star colour at the star's F606W (median, 16-84%): HIP 46120 {cal.dc_46120.median():+.4f} "
      f"({cal.dc_46120.quantile(.16):+.4f} to {cal.dc_46120.quantile(.84):+.4f}); HIP 106924 {cal.dc_106924.median():+.4f} "
      f"({cal.dc_106924.quantile(.16):+.4f} to {cal.dc_106924.quantile(.84):+.4f})")
print(f"effective number of runs, sum(p)^2 / sum(p^2): {cal.p.sum() ** 2 / (cal.p ** 2).sum():.1f}")

# which Monte Carlo inputs set the calibration-star colour offset
par = pd.read_csv(f"{B}/m15_final_refined/posterior_scaled/posterior_weights.csv").drop_duplicates("simulation_id")
pcols = ["FeH", "Y", "alphafe", "CMIXLA", "FGRY", "FGRZ", "KTTAU", "ALPHAE", "alexcoef", "opalcoef2"]
j = cal.merge(par[["simulation_id"] + pcols], on="simulation_id")
print("Spearman correlation of the HIP 46120 colour offset with each input: "
      + ", ".join(f"{c} {j.dc_46120.corr(j[c], method='spearman'):+.2f}" for c in pcols))


def wq(x, w, q):
    o = np.argsort(x)
    c = np.cumsum(w[o]) / w.sum()
    return np.interp(q, c, x[o])


print("\nposterior age (Gyr): mean, median, 16-84%; before -> after calibration weights; effective runs after")
for d in ["m15_final_refined/posterior_strict", "m15_final_refined/posterior_shifted", "m15_final_refined/posterior_scaled",
          "m15_final_refined/posterior_shifted_DB", "m15_final_refined/posterior_scaled_DB",
          "m15_yingbox_nozp_refined/posterior_shifted", "m15_yingbox_nozp_refined/posterior_scaled"]:
    w = pd.read_csv(f"{B}/{d}/posterior_weights.csv", usecols=["simulation_id", "age_gyr", "P_KS_normalized"])
    w = w.merge(cal[["simulation_id", "p"]], on="simulation_id")
    a = w.age_gyr.values
    out = []
    for ww in (w.P_KS_normalized.values, w.P_KS_normalized.values * w.p.values):
        out.append(f"{np.average(a, weights=ww):.2f}, {wq(a, ww, .5):.2f}, {wq(a, ww, .16):.2f}-{wq(a, ww, .84):.2f}")
    r = (w.P_KS_normalized * w.p).groupby(w.simulation_id).sum()
    print(f"  {d.replace('_refined/posterior', ''):34s} {out[0]}  ->  {out[1]}   (runs {r.sum() ** 2 / (r ** 2).sum():.0f})")

print("\nbest KS at each age, all runs vs runs with p > 0.05 (and their implied mu0 = mu - 2.82 E)")
for base in ["m15_final_refined", "m15_yingbox_nozp_refined"]:
    g = pd.read_csv(f"{B}/{base}/combined_2d_ks_results.csv",
                    usecols=["simulation_id", "age_gyr", "distance_modulus", "extinction_ebv", "ks_statistic"])
    g = g[g.age_gyr < g.simulation_id.map(lim) - 1e-6]
    good = g[g.simulation_id.isin(cal.simulation_id[cal.p > 0.05])]
    print(f"  {base}: {good.simulation_id.nunique()} runs with p > 0.05")
    ca = g.groupby("age_gyr").ks_statistic.min()
    gb = good.loc[good.groupby("age_gyr").ks_statistic.idxmin()].set_index("age_gyr")
    for age in np.arange(12.0, 15.01, 0.4):
        k = np.isclose(gb.index, age)
        if k.any():
            r = gb[k].iloc[0]
            print(f"    {age:4.1f}: all {ca[np.isclose(ca.index, age)].iloc[0]:.4f} | p>0.05 {r.ks_statistic:.4f} "
                  f"(run {int(r.simulation_id)}, mu {r.distance_modulus:.2f}, E {r.extinction_ebv:.3f}, "
                  f"mu0 {r.distance_modulus - 2.8235 * r.extinction_ebv:.3f})")
    print(f"    best: all {ca.idxmin():.1f} Gyr ({ca.min():.4f}); p>0.05 {gb.ks_statistic.idxmin():.1f} Gyr ({gb.ks_statistic.min():.4f})")
