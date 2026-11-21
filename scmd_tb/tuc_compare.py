#!/usr/bin/env python3
"""Observed 47 Tuc CMD against sCMDs, under the current preset's selection.

    tuc_compare.py OUT.png [--errors pairs] [--binary F] RUN:AGE:MU:EVI [RUN:AGE:MU:EVI ...]

Each sCMD is regenerated from its seed (the catalogue the grid scores for that
cell) and shown next to the observed fit stars, with a Hess residual map
(observed - sCMD) / sqrt(observed + sCMD) in 0.02 x 0.05 mag cells, and its
2D KS and JS from stage 3.
"""
import sys
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, "/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox")
from isochronetoolbox.scmd import context, preset, regenerate          # noqa: E402

ISO = "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso/{run}/isos/isochrones_example.txt"
ERRORS = "independent"
BINARY = 0.0
_ctx = None


def one(cell):
    global _ctx
    if _ctx is None:
        _ctx = context(preset("47tuc", iso=ISO, ast_errors=ERRORS, binary_fraction=BINARY))
    run, age, mu, e = cell
    cat = regenerate(_ctx, run, age, mu, e)
    return cell, _ctx.points(cat), _ctx.stage3(cat)


def main():
    global ERRORS, BINARY
    out, args = sys.argv[1], sys.argv[2:]
    while args[0].startswith("--"):
        if args[0] == "--errors":
            ERRORS = args[1]
        elif args[0] == "--binary":
            BINARY = float(args[1])
        args = args[2:]
    cells = [(int(r), float(a), float(m), float(e)) for r, a, m, e in (c.split(":") for c in args)]
    with Pool(min(len(cells), 6)) as p:
        res = p.map(one, cells)
    ctx = context(preset("47tuc", iso=ISO))
    obs = ctx.obs
    xb, yb = np.arange(0.35, 0.951, 0.02), np.arange(15.0, 19.051, 0.05)
    ho, _, _ = np.histogram2d(obs[:, 0], obs[:, 1], [xb, yb])
    n = len(res)
    fig, ax = plt.subplots(2, n + 1, figsize=(3.6 * (n + 1), 10), sharex=True, sharey=True)
    ax[0, 0].scatter(obs[:, 0], obs[:, 1], s=0.15, c="0.25", alpha=0.35, linewidths=0, rasterized=True)
    ax[0, 0].set_title(f"Observed 47 Tuc\n{len(obs)} fit stars", fontsize=10)
    ax[1, 0].pcolormesh(xb, yb, np.log10(ho.T + 1), cmap="Greys", shading="flat")
    ax[1, 0].set_title("observed Hess (log N)", fontsize=10)
    for k, ((run, age, mu, e), pts, st) in enumerate(res, 1):
        ax[0, k].scatter(pts[:, 0], pts[:, 1], s=0.15, c="#c0392b", alpha=0.35, linewidths=0, rasterized=True)
        ax[0, k].set_title(f"run {run}  {age:.1f} Gyr\nmu {mu:.2f}  E(V-I) {e:.3f}\n"
                           f"KS {st['ks_statistic']:.4f}  JS {st['js_divergence']:.4f}", fontsize=9)
        hs, _, _ = np.histogram2d(pts[:, 0], pts[:, 1], [xb, yb])
        tot = ho + hs
        z = np.where(tot > 0, (ho - hs) / np.sqrt(np.maximum(tot, 1)), np.nan)
        im = ax[1, k].pcolormesh(xb, yb, z.T, cmap="RdBu_r", vmin=-6, vmax=6, shading="flat")
        ax[1, k].set_title("(obs - sCMD) / sqrt(obs + sCMD)", fontsize=9)
    for a in ax.ravel():
        a.set_xlim(0.35, 0.95)
        a.set_ylim(19.05, 14.95)
    for a in ax[1]:
        a.set_xlabel("F606W - F814W")
    ax[0, 0].set_ylabel("F606W")
    ax[1, 0].set_ylabel("F606W")
    fig.colorbar(im, ax=ax[1, :].tolist(), shrink=0.8, pad=0.01, label="red: too few sCMD stars   blue: too many")
    fig.suptitle("47 Tuc: observed fit stars against sCMDs, both selected with the ridge cut "
                 f"(own ridge per CMD, 0.15 normalised); AST errors: {ERRORS}; binary fraction {BINARY}", fontsize=12)
    fig.savefig(out, dpi=110, bbox_inches="tight")
    for (cell, _, st) in res:
        print(cell, f"KS {st['ks_statistic']:.4f} JS {st['js_divergence']:.4f} n {st['n_simulated']}")


if __name__ == "__main__":
    main()
