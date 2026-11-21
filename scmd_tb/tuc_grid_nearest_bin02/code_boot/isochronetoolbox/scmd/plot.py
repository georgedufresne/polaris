"""Observed CMD next to a synthetic one."""

from __future__ import annotations

import numpy as np


def cmd_pair(obs: np.ndarray, sim: np.ndarray, title: str, path: str,
             window=((0.9, 1.7), (16.25, 18.0)), labels=("BP - RP", "G"),
             xlim=(0.85, 1.75), ylim=(18.6, 15.4)) -> None:
    """Observed CMD (left) and an sCMD (right) on identical axes; both are
    (colour, magnitude) arrays.  The dashed box is the 2D KS window."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    (c0, c1), (m0, m1) = window
    big = max(len(obs), len(sim)) > 5000
    fig, axes = plt.subplots(1, 2, figsize=(11, 6.2), sharey=True)
    for ax, (pts, col, name) in zip(axes, [(obs, "black", "Observed"), (sim, "#c0392b", "Synthetic")]):
        x, y = pts[:, 0], pts[:, 1]
        inside = (x >= c0) & (x <= c1) & (y >= m0) & (y <= m1)
        ax.scatter(x, y, s=0.3 if big else 2, color=col, alpha=0.35 if big else 0.6, linewidths=0, rasterized=True)
        ax.plot([c0, c1, c1, c0, c0], [m0, m0, m1, m1, m0], ls="--", lw=0.8, color="0.45")
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_xlabel(labels[0])
        ax.set_title(f"{name} CMD  (N = {len(x)}, {inside.sum()} in KS window)", fontsize=11)
        ax.grid(True, ls=":", alpha=0.4)
    axes[0].set_ylabel(labels[1])
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
