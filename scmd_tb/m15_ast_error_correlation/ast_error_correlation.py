#!/usr/bin/env python3
"""Are the F606W and F814W errors of M15 stars independent? Two checks on M. Ying's M15 data.

1. Artificial stars (Files_from_Martin/M15_data/M15artstars.dat, 100,000 stars, instrumental magnitudes): for recovered
   stars (within 0.5 px and 0.75 mag in both bands, as MakeTestCMD and our pipeline), the correlation of the two bands'
   errors (input - output), in bins of F606W (apparent = instrumental + 31.6176) and distance from the centre. Robust
   numbers: Spearman rho, and Pearson r of values winsorised at 5 x 1.4826 MAD.
2. Main-sequence colour width, F606W 19.4-20.6 (below the turnoff at 19.05), in the raw catalogue (M15.dat, no
   fiducial cut, so the distribution is not truncated): residuals in F606W-F814W from a smooth ridge (running median, then
   a degree-5 polynomial), core width 1.4826 x MAD per 0.2 mag. Predicted widths: each observed star gets an artificial
   star from its own (F606W, radius) cell, and its colour error is either that star's (dV - dI) ("joint", correlation
   kept) or dV from one artificial star and dI from another in the same cell ("independent", as MakeTestCMD draws
   them). The observed width^2 = intrinsic^2 + error^2, so an error model that predicts more width than observed is
   ruled out.
Output: logs/ast_error_correlation.out, ms_width.csv."""
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

M = "/dartfs-hpc/scratch/f007058/Files_from_Martin/M15_data"
ZPV, ZPI = 31.6176, 30.8912


def centre(x, y, guess=(3000.0, 3000.0), radius=200.0):
    near = (x - guess[0]) ** 2 + (y - guess[1]) ** 2 <= radius ** 2
    return x[near].mean(), y[near].mean()


def mad(x):
    return 1.4826 * np.median(np.abs(x - np.median(x)))


def wins_r(a, b):
    def w(x):
        m, s = np.median(x), mad(x)
        return np.clip(x, m - 5 * s, m + 5 * s)
    return np.corrcoef(w(a), w(b))[0, 1]


a = pd.read_csv(f"{M}/M15artstars.dat", sep=r"\s+", comment="#", header=None,
                names=["xi", "yi", "vi", "ii", "xo", "yo", "vo", "io"])
rec = (np.hypot(a.xi - a.xo, a.yi - a.yo) < 0.5) & ((a.vi - a.vo).abs() < 0.75) & ((a.ii - a.io).abs() < 0.75)
a = a[rec].copy()
a["dv"], a["di"] = a.vi - a.vo, a.ii - a.io
a["V"] = a.vi + ZPV
obs = pd.read_csv(f"{M}/M15.dat", sep=r"\s+", usecols=["x", "y", "Vvega", "VIvega"])
xc, yc = centre(obs.x.values, obs.y.values)
a["r"] = np.hypot(a.xi - xc, a.yi - yc)
obs["r"] = np.hypot(obs.x - xc, obs.y - yc)
print(f"artificial stars recovered: {len(a)} of 100000; centre ({xc:.1f}, {yc:.1f}) px")
rq = np.quantile(a.r, [0, 1 / 3, 2 / 3, 1])
print("\n1. correlation of the F606W and F814W errors (input - output) of the same artificial star")
print("   F606W      radius px        n    sd(dV)   sd(dI)  Spearman  Pearson(w)  colour width joint / independent")
for v0 in np.arange(17.0, 21.5, 0.5):
    for k in range(3):
        s = a[(a.V >= v0) & (a.V < v0 + 0.5) & (a.r >= rq[k]) & (a.r < rq[k + 1])]
        if len(s) < 200:
            continue
        rng = np.random.default_rng(0)
        ind = s.dv.values - rng.permutation(s.di.values)
        print(f"   {v0:4.1f}-{v0 + .5:4.1f}  {rq[k]:5.0f}-{rq[k + 1]:5.0f}  {len(s):6d}  {mad(s.dv):.4f}  {mad(s.di):.4f}  "
              f"{spearmanr(s.dv, s.di)[0]:+.2f}     {wins_r(s.dv.values, s.di.values):+.2f}      "
              f"{mad(s.dv - s.di):.4f} / {mad(ind):.4f}")
s = a[(a.V > 16.6) & (a.V < 20.7)]
print(f"   all, F606W 16.6-20.7: Spearman {spearmanr(s.dv, s.di)[0]:+.2f}, Pearson(w) {wins_r(s.dv.values, s.di.values):+.2f}")

# 2. main-sequence colour width
ms = obs[(obs.Vvega > 19.3) & (obs.Vvega < 20.7) & (obs.VIvega > 0.3) & (obs.VIvega < 1.0)].copy()
for _ in range(3):                                           # ridge: running median, then polynomial, clipped
    b = np.arange(19.3, 20.71, 0.05)
    mid = 0.5 * (b[1:] + b[:-1])
    med = [np.median(ms.VIvega[(ms.Vvega >= lo) & (ms.Vvega < hi)]) for lo, hi in zip(b[:-1], b[1:])]
    p = np.polyfit(mid, med, 5)
    ms["res"] = ms.VIvega - np.polyval(p, ms.Vvega)
    ms = ms[ms.res.abs() < 0.15]
cells_v = np.arange(19.4, 20.61, 0.2)
rq2 = np.quantile(a.r, np.linspace(0, 1, 11))
rng = np.random.default_rng(1)
print("\n2. main-sequence colour width (core, 1.4826 x MAD of the residual from the ridge)")
print("   F606W       n_obs   observed   predicted: joint  independent   implied intrinsic (joint)")
rows = []
for lo in cells_v[:-1]:
    o = ms[(ms.Vvega >= lo) & (ms.Vvega < lo + 0.2)]
    cj, ci = [], []
    for k in range(10):
        oo = o[(o.r >= rq2[k]) & (o.r < rq2[k + 1])]
        s = a[(a.V >= lo) & (a.V < lo + 0.2) & (a.r >= rq2[k]) & (a.r < rq2[k + 1])]
        if len(oo) == 0 or len(s) < 20:
            continue
        j = rng.integers(0, len(s), len(oo))
        j2 = rng.integers(0, len(s), len(oo))
        cj.append(s.dv.values[j] - s.di.values[j])
        ci.append(s.dv.values[j] - s.di.values[j2])
    cj, ci = np.concatenate(cj), np.concatenate(ci)
    w_o, w_j, w_i = mad(o.res), mad(cj), mad(ci)
    intr = np.sqrt(max(w_o ** 2 - w_j ** 2, 0))
    flag = "  <- independent predicts more than observed" if w_i > w_o else ""
    print(f"   {lo:4.1f}-{lo + .2:4.1f}  {len(o):6d}   {w_o:.4f}      {w_j:.4f}      {w_i:.4f}        {intr:.4f}{flag}")
    rows.append(dict(f606w_lo=lo, f606w_hi=lo + 0.2, n_obs=len(o), observed=w_o, joint=w_j, independent=w_i, intrinsic_joint=intr))
pd.DataFrame(rows).to_csv("/dartfs-hpc/scratch/f007058/scmd_tb/m15_ast_error_correlation/ms_width.csv", index=False)
