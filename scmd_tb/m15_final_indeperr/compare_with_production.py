#!/usr/bin/env python3
"""Production (../m15_final_refined; one artificial star's errors in both filters) against this run
(../m15_final_indeperr_refined; independent per-filter errors, the advisor's recommendation), otherwise identical:
best KS by age in sigma above each run's own bootstrap null, the overall best fit, and the strict / shifted / scaled
posterior age, DM and E. Output: logs/compare_with_production.out (via run.sh)."""
import json
import os

import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
LIM = {int(k): v for k, v in json.load(open("/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json")).items()}
RUNS = [("production (joint errors)", f"{B}/m15_final_refined"),
        ("independent per-filter errors", f"{B}/m15_final_indeperr_refined")]

curves = {}
for name, d in RUNS:
    g = pd.read_csv(f"{d}/combined_2d_ks_results.csv",
                    usecols=["simulation_id", "age_gyr", "distance_modulus", "extinction_ebv", "ks_statistic"])
    g = g[g.age_gyr < g.simulation_id.map(LIM) - 1e-6]
    n = json.load(open(f"{d}/bootstrap_summary.json"))
    m, s = n["ks_mean"], n["ks_std"]
    b = g.loc[g.ks_statistic.idxmin()]
    print(f"{name}: null {m:.5f} +/- {s:.5f}; best KS {b.ks_statistic:.5f} ({(b.ks_statistic - m) / s:.1f} sigma) at run "
          f"{int(b.simulation_id)}, {b.age_gyr:.1f} Gyr, mu {b.distance_modulus:.3f}, E {b.extinction_ebv:.3f}; "
          f"cells within 3 sigma: {(g.ks_statistic <= m + 3 * s).sum()}")
    curves[name] = (g.groupby("age_gyr").ks_statistic.min() - m) / s
    for v in ("strict", "shifted", "scaled"):
        p = f"{d}/posterior_{v}/posterior_summary.csv"
        if not os.path.exists(p):
            print(f"  {v:8s} no posterior")
            continue
        t = pd.read_csv(p).set_index("parameter")
        print(f"  {v:8s} age {t.loc['age', 'mean']:.2f} (median {t.loc['age', 'median']:.2f}, 16-84% "
              f"{t.loc['age', 'p16']:.2f}-{t.loc['age', 'p84']:.2f}); DM {t.loc['dm', 'median']:.3f}, "
              f"E {t.loc['ebv', 'median']:.3f}; {int(t.loc['age', 'runs'])} runs")
c = pd.DataFrame(curves)
print("\nbest fit at each age, sigma above own null:")
print(c[(c.index >= 11.0) & (c.index <= 15.0) & (((c.index * 10).round() % 4) == 0)].round(1).to_string())
