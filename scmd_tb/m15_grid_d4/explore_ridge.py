#!/usr/bin/env python3
"""M15 observed CMD next to 47 Tuc's: turnoff colour, sequence colours at the window edges,
and the ridge cut with the 47 Tuc settings (box shifted by the turnoff-colour difference)."""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "code")
from isochronetoolbox.scmd.ridge import RidgeCut, Tracer, turnoff_color

M = "/dartfs-hpc/scratch/f007058/Files_from_Martin"
WARP = (0.04, ((0.055, -0.05), (0.095, -0.15), (0.135, -0.30), (0.19, -0.45)), -0.5)


def load(p):
    r = pd.read_csv(p, sep=r"\s+")
    return r.VIvega.values, r.Vvega.values


def seq_colour(c, v, m, box):
    s = (np.abs(v - m) < 0.1) & (c > box[0]) & (c < box[1])
    return np.median(c[s]), s.sum()


fig, axes = plt.subplots(1, 2, figsize=(13, 8))
for ax, (name, path, sgb, box, excl) in zip(axes, (
        ("47 Tuc", f"{M}/47Tuc_data/47Tuc.dat", 17.02, (0.35, 0.95), (((0.40, 0.60), (15.5, 16.75)),)),
        ("M15", f"{M}/M15_data/M15.dat", 18.64, None, None))):
    c, v = load(path)
    lo, hi = sgb - 2, sgb + 2
    w = (v > lo) & (v < hi) & (c > 0.0) & (c < 1.5)
    to = turnoff_color(c[w], v[w])
    if box is None:                       # M15: the 47 Tuc box and exclusion shifted by the turnoff difference
        d = to - to47
        box = (round(0.35 + d, 3), round(0.95 + d, 3))
        excl = (((round(0.40 + d, 3), round(0.60 + d, 3)), (round(15.5 + sgb - 17.02, 2), round(16.75 + sgb - 17.02, 2))),)
    else:
        to47 = to
    top, nt = seq_colour(c, v, lo + 0.05, (to, 1.5))
    bot, nb = seq_colour(c, v, hi - 0.05, (to - 0.1, 1.5))
    print(f"{name}: window F606W {lo:.2f}-{hi:.2f}, {w.sum()} stars; turnoff colour {to:.3f}; "
          f"RGB colour at the top {top:.3f}; MS colour at the bottom {bot:.3f}; box {box}; exclude {excl}")
    tracer = Tracer(warp=WARP, warp_from_turnoff=True, smoothing=5e-5)
    rc, keep = RidgeCut.from_observed(tracer, c, v, (lo, hi), box, 0.15, 0.3, None, excl)
    print(f"   scales {rc.scales[0]:.3f} x {rc.scales[1]:.3f} (colour x mag): width 0.15 = "
          f"{0.15 * rc.scales[0]:.3f} in colour, {0.15 * rc.scales[1]:.3f} in mag; kept {keep.sum()} of {rc.in_box(c, v).sum()} in the box")
    ax.scatter(c[w], v[w], s=0.2, c="0.75", rasterized=True)
    ax.scatter(c[keep], v[keep], s=0.2, c="#2a78d6", rasterized=True)
    ax.plot(rc.ridge.color, rc.ridge.mag, c="#eb6834", lw=1.2)
    for (c0, c1), (m0, m1) in excl:
        ax.add_patch(plt.Rectangle((c0, m0), c1 - c0, m1 - m0, fill=False, ec="k", ls="--"))
    ax.axvline(box[0], c="k", lw=0.6); ax.axvline(box[1], c="k", lw=0.6)
    ax.set_xlim(to - 0.35, to + 0.75); ax.set_ylim(hi + 0.3, lo - 0.3)
    ax.set_title(f"{name}: kept {keep.sum()} (blue), ridge (orange)")
    ax.set_xlabel("F606W - F814W"); ax.set_ylabel("F606W")
    np.save(f"inputs/{name.replace(' ', '')}_keep.npy", keep)
fig.tight_layout()
fig.savefig("analysis_ridge_m15_vs_47tuc.png", dpi=130)
