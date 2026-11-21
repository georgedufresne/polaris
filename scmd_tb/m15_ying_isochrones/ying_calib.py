#!/usr/bin/env python3
"""The calibration-star test (Ying et al. 2023/2025; same implementation as ../m15_calib_stars/calib_weights.py)
applied to Ying's published M15 isochrones (M15_MC_iso.h5) for our 500 runs, next to our isochrones' result
(../m15_calib_stars/calib_chi2.csv). His isochrones carry his own absolute F606W and F606W-F814W, so this tests
his isochrones with his colour transformation. Ages 8-16 Gyr (his 0.2 Gyr grid). Output: logs/ying_calib.out,
ying_calib_chi2.csv."""
import h5py
import numpy as np
import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
H = f"{B}/m15_ying_isochrones"
YF = h5py.File("/dartfs-hpc/scratch/f007058/M15/zenodo_ying2025/M15_MC_iso.h5", "r")
STARS = {"HIP 46120": (5.7867, 0.0026, 0.566, 0.002), "HIP 106924": (6.0406, 0.0037, 0.601, 0.005)}
ours = pd.read_csv(f"{B}/m15_calib_stars/calib_chi2.csv").set_index("simulation_id")


def dist2(x, y, x0, y0):
    ax, ay, bx, by = x[:-1], y[:-1], x[1:], y[1:]
    dx, dy = bx - ax, by - ay
    t = np.clip(((x0 - ax) * dx + (y0 - ay) * dy) / np.maximum(dx * dx + dy * dy, 1e-30), 0, 1)
    return np.min((ax + t * dx - x0) ** 2 + (ay + t * dy - y0) ** 2)


rows = []
for run in ours.index:
    g = YF[f"{run}/isochrones/cmd"]
    best = None
    for age in g:
        a = int(age) / 1000
        if a < 8.0 - 1e-6 or a > 16.0 + 1e-6:
            continue
        d = g[age]
        cols = [x.decode() for x in d["block0_items"][()]]
        t = pd.DataFrame(d["block0_values"][()], columns=cols).sort_values("EEP")
        v, c = t.F606W.values, t.F814W.values                 # 'F814W' column holds F606W-F814W
        chi2, dc = 0.0, {}
        for name, (M, sM, C, sC) in STARS.items():
            chi2 += dist2(c / sC, v / sM, C / sC, M / sM)
            ms = t[(v > 4.5) & (v < 7.5)].sort_values("F606W")
            dc[name] = np.interp(M, ms.F606W.values, ms.F814W.values) - C
        if best is None or chi2 < best[1]:
            best = (a, chi2, dc)
    a, chi2, dc = best
    rows.append(dict(simulation_id=run, chi2_min=chi2, age_at_min=a, p=np.exp(-chi2 / 2),
                     dc_46120=dc["HIP 46120"], dc_106924=dc["HIP 106924"]))
y = pd.DataFrame(rows).set_index("simulation_id")
y.to_csv(f"{H}/ying_calib_chi2.csv")
for name, t in (("his isochrones", y), ("our isochrones", ours.loc[y.index])):
    print(f"{name}: chi2_min median {t.chi2_min.median():.1f}; chi2 < 2: {(t.chi2_min < 2).mean():.0%}; p > 0.05: "
          f"{(t.p > 0.05).mean():.0%}; colour offset at the stars' F606W (isochrone - star), median (16-84%): "
          f"HIP 46120 {t.dc_46120.median():+.4f} ({t.dc_46120.quantile(.16):+.4f} to {t.dc_46120.quantile(.84):+.4f}), "
          f"HIP 106924 {t.dc_106924.median():+.4f}")
d = y.dc_46120 - ours.loc[y.index, "dc_46120"]
print(f"his minus our main-sequence colour at HIP 46120's F606W, same run: median {d.median():+.4f} "
      f"(16-84% {d.quantile(.16):+.4f} to {d.quantile(.84):+.4f})")
