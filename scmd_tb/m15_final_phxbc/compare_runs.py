#!/usr/bin/env python3
"""Compare the production-methodology runs, each against its own bootstrap null:
  production (our isochrones, joint errors)       ../m15_final_refined
  + PHOENIX colour tables (our models)            ../m15_final_phxbc_refined
  + independent per-filter errors                 ../m15_final_indeperr_refined
  + Ying's isochrones                          ../m15_final_yingiso_refined
  + Ying's isochrones + independent errors        ../m15_final_yingiso_indeperr_refined
and, for reference, Ying's isochrones in our pipeline without the map (his E box, joint errors), ../m15_yingiso_nozp_refined.
Prints the best fit (sigma above the null, cells within 3 sigma), strict / shifted / scaled posteriors (age, DM, E),
and the best fit at each age in sigma above the run's own null. Output: logs/compare_runs.out (via run.sh)."""
import json
import os

import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
D4 = "/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json"
YL = "/dartfs-hpc/scratch/f007058/M15/iso_ying2025/age_limits_ying.json"
RUNS = [("production", f"{B}/m15_final_refined", D4),
        ("+ PHOENIX BCs", f"{B}/m15_final_phxbc_refined", D4),
        ("+ indep. errors", f"{B}/m15_final_indeperr_refined", D4),
        ("+ Ying iso", f"{B}/m15_final_yingiso_refined", YL),
        ("+ Ying iso + indep.", f"{B}/m15_final_yingiso_indeperr_refined", YL),
        ("Ying iso, no map (ref)", f"{B}/m15_yingiso_nozp_refined", YL)]

curves = {}
for name, d, limf in RUNS:
    if not os.path.exists(f"{d}/bootstrap_summary.json"):
        print(f"{name}: not finished")
        continue
    lim = {int(k): v for k, v in json.load(open(limf)).items()}
    g = pd.read_csv(f"{d}/combined_2d_ks_results.csv",
                    usecols=["simulation_id", "age_gyr", "distance_modulus", "extinction_ebv", "ks_statistic"])
    g = g[g.age_gyr < g.simulation_id.map(lim) - 1e-6]
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
ages = [a for a in c.index if 11.0 - 1e-6 <= a <= 15.0 + 1e-6 and abs(a * 10 % 4) < 1e-6]
print("\nbest fit at each age, sigma above own null:")
print(c.loc[ages].round(1).to_string())
