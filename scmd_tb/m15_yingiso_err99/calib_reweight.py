#!/usr/bin/env python3
"""Ying's full weighting on this run: the posterior weights (P_KS_normalized) times his calibration-star weight per run,
w_cal = exp(-chi2_min/2), computed for HIS isochrones (../m15_ying_isochrones/ying_calib_chi2.csv; recipe of Ying et al.
2023/2025, see ../m15_calib_stars/). Output: logs/calib_reweight.out."""
import numpy as np
import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
cal = pd.read_csv(f"{B}/m15_ying_isochrones/ying_calib_chi2.csv")


def wq(x, w, q):
    o = np.argsort(x)
    c = np.cumsum(w[o]) / w.sum()
    return np.interp(q, c, x[o])


print("posterior age (Gyr): mean, median, 16-84%; KS weights only -> times calibration weights (his isochrones); effective runs")
for v in ("shifted", "scaled"):        # no strict posterior: no cell within 3 sigma
    w = pd.read_csv(f"{B}/m15_yingiso_err99_refined/posterior_{v}/posterior_weights.csv",
                    usecols=["simulation_id", "age_gyr", "P_KS_normalized"]).merge(cal[["simulation_id", "p"]], on="simulation_id")
    a = w.age_gyr.values
    out = [f"{np.average(a, weights=x):.2f}, {wq(a, x, .5):.2f}, {wq(a, x, .16):.2f}-{wq(a, x, .84):.2f}"
           for x in (w.P_KS_normalized.values, w.P_KS_normalized.values * w.p.values)]
    r = (w.P_KS_normalized * w.p).groupby(w.simulation_id).sum()
    print(f"  {v:8s} {out[0]}  ->  {out[1]}   (runs {r.sum() ** 2 / (r ** 2).sum():.0f}; "
          f"weight-carrying runs passing p > 0.05: {(w.drop_duplicates('simulation_id').p > 0.05).mean():.0%})")
