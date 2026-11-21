#!/usr/bin/env python3
"""Why do Ying's M15 isochrones differ from ours when the inputs are the same? Colour transformation or stellar models?

Ying et al. (2025) convert with Dartmouth PHOENIX bolometric corrections at the nearest [alpha/Fe]. Production converts
our DSEP4 models (log L, log Teff, log g) with one MIST ACS/WFC table at [Fe/H] -2.25 (BCTable, as in the pipeline).
Here the same model points are converted three ways and compared with his published isochrones:
  A  our models + MIST (production; ../m15_final/code BCTable)
  B  our models + PHOENIX (ChaboyerB/Catherine/dsep3/iso/phx/Z*.HST_ACSWF: BCv, V-F606W, V-F814W on (Teff, log g);
     linear in (Teff, log g) on the table points, linear in [Fe/H] between -2.0 and -2.5 (the run's [Fe/H]),
     nearest [alpha/Fe] of -0.2..0.8, M_bol = 4.75 - 2.5 log L as in BCTable; in these tables BCv = M_V - M_bol (the Sun,
     5800 K, log g 4.5, solar: +0.082, i.e. M_V 4.83), so M_V = M_bol + BCv; points outside the table are dropped)
  C  his isochrones (Zenodo 13998637)
1. The two transformations directly: PHOENIX minus MIST F606W-F814W at typical M15 (Teff, log g) points.
2. Same 29 runs as ../m15_ying_isochrones: his minus A and his minus B in turnoff, main-sequence (M_F606W 4.5, 5.8) and
   RGB colour, and the SGB-turnoff F606W gap, at 13.4 and 14.2 Gyr.
3. Equivalent age (as ../m15_ying_isochrones/compare_isochrones.py): which of our ages, in A and in B, best matches his
   13.4 Gyr isochrone within +/-2 mag of its SGB point, with free F606W and colour shifts.
4. Calibration stars (HIP 46120, HIP 106924; recipe of ../m15_calib_stars) for B on the same runs, next to A and C.
If B reproduces C, the difference is the colour transformation; if not, it is in the models.
Output: logs/bc_test.out, bc_test.csv."""
import glob
import json
import sys

import h5py
import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator
from scipy.optimize import minimize
from scipy.spatial import cKDTree

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
H = f"{B}/m15_isochrone_bc_test"
PHX = "/dartfs-hpc/rc/lab/C/ChaboyerB/Catherine/dsep3/iso/phx"
sys.path.insert(0, f"{B}/m15_final/code")
from isochronetoolbox.scmd.model import BCTable, read_isochrones  # noqa: E402

K = 6.0
YF = h5py.File("/dartfs-hpc/scratch/f007058/M15/zenodo_ying2025/M15_MC_iso.h5", "r")
mist = BCTable(f"{B}/m15_grid_d4/inputs/bc_mist_acs_wfc_fehm225.npz")
prev = pd.read_csv(f"{B}/m15_ying_isochrones/compare_isochrones.csv").set_index("simulation_id")
RUNS = list(prev.index)
LIM = {int(k): v for k, v in json.load(open("/dartfs-hpc/scratch/f007058/M15/dsep4/m15_age_limits_d4.json")).items()}
STARS = {"HIP 46120": (5.7867, 0.0026, 0.566, 0.002), "HIP 106924": (6.0406, 0.0037, 0.601, 0.005)}


class Phoenix:
    """Dartmouth PHOENIX BCs: (Teff, log g) -> (BCv, V-F606W, V-F814W) per ([Fe/H], [alpha/Fe]) table."""
    FEH = {"m2d0": -2.0, "m2d5": -2.5}
    AFE = {"am2": -0.2, "ap0": 0.0, "ap2": 0.2, "ap4": 0.4, "ap6": 0.6, "ap8": 0.8}

    def __init__(self):
        self.t = {}
        for fz, feh in self.FEH.items():
            for fa, afe in self.AFE.items():
                d = pd.read_csv(f"{PHX}/Z{fz}{fa}.HST_ACSWF", sep=r"\s+", comment="#", header=None).values
                self.t[(feh, afe)] = LinearNDInterpolator(d[:, :2], d[:, [3, 7, 10]])

    def absolute(self, rows, feh, afe):
        a = min(self.AFE.values(), key=lambda x: abs(x - afe))
        f = np.clip((feh - -2.0) / (-2.5 - -2.0), 0, 1)
        x = (10 ** rows["log_Teff"].values, rows["log_g"].values)
        bc = (1 - f) * self.t[(-2.0, a)](*x) + f * self.t[(-2.5, a)](*x)
        mv = 4.75 - 2.5 * rows["log_L"].values + bc[:, 0]          # BCv = M_V - M_bol in these tables
        return {"F606W": mv - bc[:, 1], "F814W": mv - bc[:, 2]}


phx = Phoenix()


def read_df(path):
    g = YF[path]
    return pd.DataFrame(g["block0_values"][()], columns=[x.decode() for x in g["block0_items"][()]])


def ying(run, age):
    d = read_df(f"{run}/isochrones/cmd/{int(round(age * 1000))}").sort_values("EEP")
    return d.F606W.values, d.F814W.values


def ours(run, conv, feh=None, afe=None):
    iso = read_isochrones(f"/dartfs-hpc/scratch/f007058/M15/iso_dsep4/{run}/isos/isochrones_example.txt")
    m = mist.absolute(iso) if conv == "A" else phx.absolute(iso, feh, afe)
    iso = iso.assign(v=m["F606W"], c=m["F606W"] - m["F814W"])
    iso = iso[np.isfinite(iso.v) & np.isfinite(iso.c)]
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
    return dict(v_to=v[i_to], c_to=c[i_to], v_sgb=vs, gap=vs - v[i_to],
                c_ms58=np.interp(5.8, v[ms][::-1], c[ms][::-1]), c_ms45=np.interp(4.5, v[ms][::-1], c[ms][::-1]),
                c_rgb=np.interp(vs - 1.5, v[rgb][::-1], c[rgb][::-1]) if len(rgb) and v[rgb].min() < vs - 1.5 else np.nan)


def dense(v, c, step=0.002):
    ok = np.isfinite(v) & np.isfinite(c)
    v, c = v[ok], c[ok]
    s = np.concatenate([[0], np.cumsum(np.hypot(K * np.diff(c), np.diff(v)))])
    t = np.arange(0, s[-1], step)
    return np.column_stack([K * np.interp(t, s, c), np.interp(t, s, v)])


def equiv_age(run, curves):
    v, c = ying(run, 13.4)
    fy = features(v, c)
    k = (v > fy["v_sgb"] - 2) & (v < fy["v_sgb"] + 2)
    target = np.column_stack([K * c[k], v[k]])
    ok = [a for a in sorted(curves) if 11.0 <= a <= 15.8 and a < LIM[run] - 1e-6]    # as compare_isochrones.py

    def fit(a):
        try:
            fo = features(*curves[a])
        except IndexError:                 # no subgiant point (isochrone cut short)
            return a, np.inf, (np.nan, np.nan)
        tree = cKDTree(dense(*curves[a]))
        f = lambda p: np.mean(tree.query(target + [K * p[1], p[0]])[0] ** 2)       # noqa: E731
        r = minimize(f, [fy["v_sgb"] - fo["v_sgb"], fy["c_to"] - fo["c_to"]], method="Nelder-Mead",
                     options=dict(xatol=1e-4, fatol=1e-10))
        return a, r.fun, r.x
    coarse = [fit(a) for a in ok[::2]]
    a0 = min(coarse, key=lambda x: x[1])[0]
    fits = coarse + [fit(a) for a in ok if abs(a - a0) < 0.25 and a not in ok[::2]]
    a, cost, x = min(fits, key=lambda x: x[1])
    return a, np.sqrt(cost), x


def dist2(x, y, x0, y0):
    ax, ay, bx, by = x[:-1], y[:-1], x[1:], y[1:]
    dx, dy = bx - ax, by - ay
    t = np.clip(((x0 - ax) * dx + (y0 - ay) * dy) / np.maximum(dx * dx + dy * dy, 1e-30), 0, 1)
    return np.min((ax + t * dx - x0) ** 2 + (ay + t * dy - y0) ** 2)


def calib(curves):
    best = np.inf
    for a, (v, c) in curves.items():
        if 8.0 - 1e-6 <= a <= 16.0 + 1e-6:
            ok = np.isfinite(v) & np.isfinite(c)
            best = min(best, sum(dist2(c[ok] / sC, v[ok] / sM, C / sC, M / sM) for M, sM, C, sC in STARS.values()))
    return best


# 1. the two transformations at typical points
print("1. PHOENIX ([Fe/H] -2.25, [a/Fe] +0.2/+0.4) minus MIST ([Fe/H] -2.25) in F606W-F814W and F606W BC:")
pts = [("lower MS", 4800, 4.8, -1.2), ("MS, M_F606W ~5.8", 5300, 4.7, -0.9), ("MS, M_F606W ~4.5", 5900, 4.6, -0.4),
       ("turnoff", 6500, 4.2, 0.4), ("subgiant", 6000, 3.8, 0.8), ("RGB base", 5300, 3.4, 1.0), ("RGB", 4900, 2.6, 1.8)]
for name, teff, logg, logl in pts:
    row = pd.DataFrame({"log_Teff": [np.log10(teff)], "log_g": [logg], "log_L": [logl]})
    m = mist.absolute(row)
    out = []
    for afe in (0.2, 0.4):
        p = phx.absolute(row, -2.25, afe)
        out.append(f"[a/Fe] {afe}: colour {(p['F606W'] - p['F814W'])[0] - (m['F606W'] - m['F814W'])[0]:+.4f}, "
                   f"F606W {p['F606W'][0] - m['F606W'][0]:+.4f}")
    print(f"   {name:18s} Teff {teff} log g {logg}: " + "; ".join(out))

# 2-4. per run
rows = []
for run in RUNS:
    var = read_df(f"{run}/var").iloc[0]
    A = ours(run, "A")
    Bc = ours(run, "B", var.FeH, var.alphafe)
    r = dict(simulation_id=run, FeH=var.FeH, alphafe=var.alphafe)
    for age in (13.4, 14.2):
        fc = features(*ying(run, age))
        for tag, cur in (("A", A), ("B", Bc)):
            f = features(*cur[age])
            for k in ("c_to", "c_ms58", "c_ms45", "c_rgb", "gap"):
                r[f"C_minus_{tag}_{k}_{age}"] = fc[k] - f[k]
    for tag, cur in (("A", A), ("B", Bc)):
        a, rms, (dv, dc) = equiv_age(run, cur)
        r.update({f"equiv_{tag}": a, f"equiv_rms_{tag}": rms, f"equiv_dc_{tag}": dc})
    r["calib_chi2_B"] = calib(Bc)
    rows.append(r)
    print(f"   run {run} [Fe/H] {var.FeH:+.2f} [a/Fe] {var.alphafe:+.2f}: his 13.4 Gyr matches A {r['equiv_A']:.1f} "
          f"(prev {prev.loc[run, 'equiv_age_for_ying_13_4']:.1f}), B {r['equiv_B']:.1f}; turnoff colour his-A "
          f"{r['C_minus_A_c_to_13.4']:+.4f}, his-B {r['C_minus_B_c_to_13.4']:+.4f}", flush=True)
d = pd.DataFrame(rows)
d.to_csv(f"{H}/bc_test.csv", index=False)

print(f"\n2. his minus ours, median (16-84%) over {len(d)} runs, at 13.4 / 14.2 Gyr:")
for k, name in (("c_to", "turnoff colour"), ("c_ms45", "MS colour at M_F606W 4.5"), ("c_ms58", "MS colour at M_F606W 5.8"),
                ("c_rgb", "RGB colour"), ("gap", "SGB - turnoff F606W gap")):
    s = []
    for tag in ("A", "B"):
        for age in (13.4, 14.2):
            x = d[f"C_minus_{tag}_{k}_{age}"]
            s.append(f"{tag} {age}: {x.median():+.4f} ({x.quantile(.16):+.4f}..{x.quantile(.84):+.4f})")
    print(f"   {name:26s} " + " | ".join(s))
print("\n3. our age that best matches his 13.4 Gyr isochrone (free shifts):")
for tag, name in (("A", "A (MIST, production)"), ("B", "B (PHOENIX)")):
    e = d[f"equiv_{tag}"]
    print(f"   {name:22s} median {e.median():.2f}, 16-84% {e.quantile(.16):.2f}-{e.quantile(.84):.2f}, range "
          f"{e.min():.1f}-{e.max():.1f}; rms of match {d[f'equiv_rms_{tag}'].median():.4f}; "
          f"|age - 13.4| <= 0.2 in {(np.abs(e - 13.4) <= 0.2001).mean():.0%} of runs")
cA = pd.read_csv(f"{B}/m15_calib_stars/calib_chi2.csv").set_index("simulation_id").loc[d.simulation_id]
cC = pd.read_csv(f"{B}/m15_ying_isochrones/ying_calib_chi2.csv").set_index("simulation_id").loc[d.simulation_id]
print(f"\n4. calibration stars on these {len(d)} runs, chi2 < 2 / median chi2: A (ours, MIST) {(cA.chi2_min < 2).mean():.0%} / "
      f"{cA.chi2_min.median():.1f}; B (ours, PHOENIX) {(d.calib_chi2_B < 2).mean():.0%} / {d.calib_chi2_B.median():.1f}; "
      f"C (his) {(cC.chi2_min < 2).mean():.0%} / {cC.chi2_min.median():.1f}")
