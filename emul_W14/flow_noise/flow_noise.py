#!/usr/bin/env python3
"""How noisy is the W v1.4 emulator, in NGC 6791 age, and is that noise what moves the W v1.4 age?

W v1.4 is a conditional neural spline flow (zuko NSF): for given physics, EEP and initial mass it is a density over
(log Teff, log g, log L, log R, central He4, He-core mass). Our isochrones use its mode (z = 0). The masses along an
isochrone come from the separate age head (deterministic MLP), so the flow's spread is a spread in POSITION at fixed mass.

For each run (NGC 6791 var file), at REF_AGE:
  A. flow spread. K realizations drawn from the flow at the reference isochrone's (EEP, mass) points:
       coherent:    one base-Gaussian draw z shared by every point (the whole curve moves together: worst case);
       independent: a fresh z at every point (per-point noise, which partly averages out in a fit).
     Each realization is fitted for age against the central (mode) isochrone grid.
  B. actual error. The DSEP3 physics-on isochrone of the same run (Isochrones_noco, what the "New" sCMD leg used) at
     8.0 / 8.3 / 8.6 Gyr, fitted the same way: the emulator's age offset with no sCMD noise.
  Self-test: central isochrones at 8.0 / 8.3 / 8.6 fitted against the grid must return their own age.

Fit: pseudo-CMD x = -5 log Teff, y = -2.5 log L (local stand-ins for colour and magnitude). Box like the sCMD 2D-KS
(G 16.25-18.0 around a turnoff at G ~17.3): from 1.05 mag brighter to 0.70 mag fainter than the data's turnoff. Each point
weighted by its share of stars (Salpeter IMF x mass interval). Free 2D shift (distance + reddening), found by
translation-only ICP on the nearest point of the densified model curve; age = parabola through the rms minimum.

  python flow_noise.py RUN [RUN ...] [--k 100] [--out output/]
Writes output/run_<RUN>.json. summarize.py combines them.
"""
import argparse
import json
import os
import sys
import time

import numpy as np
import torch
from scipy.spatial import cKDTree

S = "/dartfs-hpc/scratch/f007058/emul_W14"
sys.path.insert(0, S)
sys.path.insert(0, f"{S}/CONF1DENCE-W_fcfac3c")
from gen_W14_iso import age_table, invert, phys_for       # noqa: E402  (the exact recipe behind Isochrones_W14)
from isochrone_W import central_mode                       # noqa: E402
from predict_W import WEmulator, AgeHead                   # noqa: E402

DSEP3 = "/dartfs-hpc/scratch/f007058/Isochrones_noco/{run}/isos/isochrones_example.txt"
REF_AGE = 8.3
GRID = np.round(np.arange(7.0, 9.9001, 0.05), 3)           # central isochrone ages the fits choose from
BRIGHT, FAINT = 1.05, 0.70                                 # box around the turnoff, in mag
MASS_LO, MASS_HI, NMASS = 0.10, 2.5, 1200


def pseudo(logT, logL):
    return -5.0 * np.asarray(logT), -2.5 * np.asarray(logL)


def turnoff(x, y, logL):
    """Bluest point before the subgiant branch (log L < 1)."""
    k = np.nonzero(np.asarray(logL) < 1.0)[0]
    i = k[np.argmin(x[k])]
    return x[i], y[i]


def imf_weights(m):
    """Share of stars at each point: Salpeter dN/dm times the mass interval around the point."""
    m = np.asarray(m, float)
    dm = np.abs(np.gradient(m))
    return m ** -2.35 * np.maximum(dm, 1e-12)


def densify(x, y, step=0.001):
    s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))])
    n = max(int(s[-1] / step), 2)
    t = np.linspace(0, s[-1], n)
    return np.interp(t, s, x), np.interp(t, s, y)


class Model:
    """One central isochrone, restricted to its own turnoff box (+0.3 mag margin), densified, with a KD tree."""
    def __init__(self, age, x, y, logL):
        self.age = age
        xt, yt = turnoff(x, y, logL)
        k = (y > yt - BRIGHT - 0.3) & (y < yt + FAINT + 0.3)
        # keep the box section as one contiguous piece of the curve (MS -> TO -> SGB -> lower RGB)
        idx = np.nonzero(k)[0]
        idx = idx[(idx >= idx.min()) & (idx <= idx.max())]
        xd, yd = densify(x[idx.min():idx.max() + 1], y[idx.min():idx.max() + 1])
        self.pts = np.column_stack([xd, yd])
        self.tree = cKDTree(self.pts)
        self.to = np.array([xt, yt])


def fit_rms(model, D, w, to_data, iters=20):
    """Translation-only ICP: best 2D shift of the data onto the model curve; returns weighted rms distance (mag)."""
    d = model.to - to_data
    for _ in range(iters):
        P = D + d
        _, j = model.tree.query(P)
        r = model.pts[j] - P
        step = (w[:, None] * r).sum(0) / w.sum()
        d = d + step
        if np.hypot(*step) < 1e-6:
            break
    P = D + d
    dist, _ = model.tree.query(P)
    return float(np.sqrt((w * dist ** 2).sum() / w.sum())), d


def fit_age(models, x, y, logL, m, box=None, ito=None):
    """Age of a data curve (x, y) against the central grid: rms minimum, refined by a parabola.
    box / ito: for draws at the reference isochrone's own points, its box mask and turnoff index (a turnoff searched on a
    jittered curve picks noise spikes)."""
    if box is None:
        xt, yt = turnoff(x, y, logL)
        k = (y > yt - BRIGHT) & (y < yt + FAINT)
    else:
        xt, yt = x[ito], y[ito]
        k = box
    D = np.column_stack([x[k], y[k]]); w = imf_weights(m)[k]
    rms = np.array([fit_rms(M, D, w, np.array([xt, yt]))[0] for M in models])
    i = int(np.argmin(rms))
    a = GRID
    if 0 < i < len(a) - 1:
        y0, y1, y2 = rms[i - 1], rms[i], rms[i + 1]
        den = y0 - 2 * y1 + y2
        off = 0.5 * (y0 - y2) / den if den > 0 else 0.0
        age = a[i] + np.clip(off, -1, 1) * (a[1] - a[0])
    else:
        age = a[i]
    return float(age), float(rms[i]), bool(0 < i < len(a) - 1)


def dsep3_block(run, age):
    M = np.loadtxt(DSEP3.format(run=run), comments="#")
    k = np.abs(M[:, 1] / 1e9 - age) < 1e-3
    B = M[k]
    return B[:, 10], B[:, 7], B[:, 2]                      # log Teff, log L, initial mass


@torch.no_grad()
def draws(W, eep, mass, phys, K, coherent, rng):
    """K realizations of (log Teff, log L) at the given points, from the flow's own density."""
    ctx = W._context(phys, eep, mass)
    T = W.flow(ctx).transform
    nf = W.cfg["n_features"]; n = len(eep)
    out = []
    for _ in range(K):
        z = rng.standard_normal(nf) * np.ones((n, 1)) if coherent else rng.standard_normal((n, nf))
        s = T.inv(torch.as_tensor(z, dtype=torch.float32)).cpu().numpy()
        o = W._physical(s)
        out.append((o["log_Teff"], o["log_L"]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+", type=int)
    ap.add_argument("--k", type=int, default=100)
    ap.add_argument("--out", default=f"{S}/flow_noise/output")
    ap.add_argument("--seed", type=int, default=20261122)
    a = ap.parse_args()
    torch.set_num_threads(int(os.environ.get("SLURM_CPUS_PER_TASK", "4")))
    W = WEmulator(); ah = AgeHead()
    lh = W.fmt["label_header"]; default = W.fmt["label_default"]
    for run in a.runs:
        out = f"{a.out}/run_{run}.json"
        if os.path.exists(out):
            print(f"{run}: exists", flush=True); continue
        t0 = time.time(); rng = np.random.default_rng(a.seed + run)
        phys = phys_for(run, f"{S}/var", lh, default, "var")
        p22 = np.asarray(phys[:22], float)[:len(ah.theta)]
        eg, mg, LA = age_table(ah, p22, MASS_LO, MASS_HI, NMASS)
        cen = {}
        for age in GRID:
            m, e = invert(ah, p22, eg, mg, LA, float(age), MASS_LO, MASS_HI)
            o = central_mode(W, e, m, phys)
            cen[float(age)] = (e, m, o["log_Teff"], o["log_L"])
        models = []
        for age in GRID:
            e, m, lT, lL = cen[float(age)]; x, y = pseudo(lT, lL)
            models.append(Model(float(age), x, y, lL))
        res = {"run": run, "ref_age": REF_AGE, "k": a.k, "grid": [float(GRID[0]), float(GRID[-1]), 0.05]}
        # self-test
        res["self"] = {}
        for age in (8.0, 8.3, 8.6):
            e, m, lT, lL = cen[age]; x, y = pseudo(lT, lL)
            res["self"][str(age)] = fit_age(models, x, y, lL, m)
        # A. flow spread at the reference isochrone's points
        e, m, lT0, lL0 = cen[REF_AGE]
        x0, y0 = pseudo(lT0, lL0); xt, yt = turnoff(x0, y0, lL0)
        box = (y0 > yt - BRIGHT) & (y0 < yt + FAINT); ito = int(np.argmin(np.where(lL0 < 1.0, x0, np.inf)))
        for kind, coh in (("coherent", True), ("independent", False)):
            D = draws(W, e, m, phys, a.k, coh, rng)
            ages, rmss, inside = [], [], []
            for lT, lL in D:
                x, y = pseudo(lT, lL)
                ag, r, ok = fit_age(models, x, y, lL, m, box=box, ito=ito)
                ages.append(ag); rmss.append(r); inside.append(ok)
            LT = np.array([d[0] for d in D]); LL = np.array([d[1] for d in D])
            res[kind] = {"ages": ages, "rms": rmss, "inside_grid": inside,
                         "sd_logTeff_TO": float(LT[:, ito].std()), "sd_logL_TO": float(LL[:, ito].std()),
                         "median_sd_logTeff_box": float(np.median(LT[:, box].std(0))),
                         "median_sd_logL_box": float(np.median(LL[:, box].std(0))),
                         "mean_minus_mode_logTeff_TO": float(LT[:, ito].mean() - lT0[ito]),
                         "mean_minus_mode_logL_TO": float(LL[:, ito].mean() - lL0[ito])}
        # B. DSEP3 isochrone of the same run
        res["dsep3"] = {}
        for age in (8.0, 8.3, 8.6):
            lT, lL, mm = dsep3_block(run, age); x, y = pseudo(lT, lL)
            xd, yd = turnoff(x, y, lL); e0, m0, lTw, lLw = cen[age]; xw, yw = pseudo(lTw, lLw); xtw, ytw = turnoff(xw, yw, lLw)
            res["dsep3"][str(age)] = list(fit_age(models, x, y, lL, mm)) + [
                float(-(xd - xtw) / 5.0), float(-(yd - ytw) / 2.5)]     # + TO offsets DSEP3 - W (dex log Teff, log L)
        res["seconds"] = round(time.time() - t0, 1)
        tmp = out + ".tmp"
        json.dump(res, open(tmp, "w"), indent=1)
        os.replace(tmp, out)
        c = np.array(res["coherent"]["ages"]); i = np.array(res["independent"]["ages"])
        print(f"{run}: self {[round(v[0], 3) for v in res['self'].values()]}  coherent {c.mean():.3f} +/- {c.std():.3f}  "
              f"independent {i.mean():.3f} +/- {i.std():.3f}  DSEP3@8.3 -> {res['dsep3']['8.3'][0]:.3f}  "
              f"({res['seconds']:.0f} s)", flush=True)


if __name__ == "__main__":
    main()
