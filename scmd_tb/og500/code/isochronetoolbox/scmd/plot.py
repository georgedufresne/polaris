"""Observed CMD next to a synthetic one."""

from __future__ import annotations

import numpy as np
import pandas as pd


def cmd_pair(obs: pd.DataFrame, cat: np.ndarray, title: str, path: str,
             ks_window=((0.9, 1.7), (16.25, 18.0)), xlim=(0.85, 1.75), ylim=(18.6, 15.4)) -> None:
    """Observed CMD (left) and an sCMD (right) on identical axes.

    The dashed box is the window the 2D KS statistic is computed in.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    (c0, c1), (m0, m1) = ks_window
    color = obs["BP_RP"].values if "BP_RP" in obs else (obs["BPmag"] - obs["RPmag"]).values
    panels = [(color, obs["Gmag"].values, "black", "Observed"),
              (cat[:, 1] - cat[:, 2], cat[:, 0], "#c0392b", "Synthetic")]
    fig, axes = plt.subplots(1, 2, figsize=(11, 6.2), sharey=True)
    for ax, (x, y, col, name) in zip(axes, panels):
        inside = (x >= c0) & (x <= c1) & (y >= m0) & (y <= m1)
        ax.scatter(x, y, s=2, color=col, alpha=0.6, linewidths=0)
        ax.plot([c0, c1, c1, c0, c0], [m0, m0, m1, m1, m0], ls="--", lw=0.8, color="0.45")
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_xlabel("BP - RP")
        ax.set_title(f"{name} CMD  (N = {len(x)}, {inside.sum()} in KS window)", fontsize=11)
        ax.grid(True, ls=":", alpha=0.4)
    axes[0].set_ylabel("G")
    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
