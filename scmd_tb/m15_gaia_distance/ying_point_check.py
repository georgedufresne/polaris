#!/usr/bin/env python3
"""Implied true distance modulus of the best fit at each age, and the best age within Ying et al. (2025)'s own
M15 fit (Table 7: DM 15.49 +/- 0.03, E 0.08 +/- 0.01), from the existing grid + refinement cells of production
(m15_final_refined) and of Ying's box without the zero-point map (m15_yingbox_nozp_refined).
Output: logs/ying_point_check.out."""
import json

import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
lim = {int(k): v for k, v in json.load(open("/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json")).items()}
A = 2.88 / 1.02          # A_F606W per unit E(F606W-F814W), ACS
print("Ying Table 7: M15 DM 15.49+-0.03, E 0.08+-0.01, 13.39+-0.39 Gyr -> implied mu0", round(15.49 - A * 0.08, 3),
      "(Gaia 15.149+-0.020)")
for base in ["m15_final_refined", "m15_yingbox_nozp_refined"]:
    g = pd.read_csv(f"{B}/{base}/combined_2d_ks_results.csv",
                    usecols=["simulation_id", "age_gyr", "distance_modulus", "extinction_ebv", "ks_statistic"])
    g = g[g.age_gyr < g.simulation_id.map(lim) - 1e-6]
    g["mu0"] = g.distance_modulus - A * g.extinction_ebv
    print(f"\n{base}")
    b = g.loc[g.groupby("age_gyr").ks_statistic.idxmin()].set_index("age_gyr")
    b = b[(b.index >= 12.0) & (b.index <= 15.0)]
    print("  best fit per age: age, KS, mu, E, implied mu0")
    for a, r in b.iloc[::2].iterrows():
        print(f"   {a:4.1f}  {r.ks_statistic:.4f}  {r.distance_modulus:.2f}  {r.extinction_ebv:.3f}  {r.mu0:.3f}")
    for name, lo_mu, hi_mu, lo_e, hi_e in (("Ying's point +-1 sigma (mu 15.46-15.52, E 0.07-0.09)", 15.455, 15.525, 0.0695, 0.0905),
                                           ("Ying's point +-2 sigma (mu 15.43-15.55, E 0.06-0.10)", 15.425, 15.555, 0.0595, 0.1005)):
        s = g[g.distance_modulus.between(lo_mu, hi_mu) & g.extinction_ebv.between(lo_e, hi_e)]
        c = s.groupby("age_gyr").ks_statistic.min()
        c = c[(c.index >= 12.0) & (c.index <= 15.0)]
        print(f"  {name}: {len(s)} cells; best {c.idxmin():.1f} Gyr KS {c.min():.4f}")
        print("   ", "  ".join(f"{a:.1f}:{v:.4f}" for a, v in c.iloc[::2].items()))
