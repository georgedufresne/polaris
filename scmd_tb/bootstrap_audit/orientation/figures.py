"""Figures for the orientation audit (static PNGs from the saved per-point arrays and the summary table)."""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

H = "/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/orientation"
INK, MUTED, GRID = "#1f1f1e", "#6b6b66", "#e4e3df"
DIV = LinearSegmentedColormap.from_list("div", ["#104281", "#3987e5", "#f0efec", "#e66767", "#9e2323"])
NAMES = {"bb": "bluer & brighter (pipeline)", "bf": "bluer & fainter (M. Ying)",
         "rb": "redder & brighter", "rf": "redder & fainter"}
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.titlesize": 9, "figure.facecolor": "#fcfcfb",
                     "axes.facecolor": "#fcfcfb"})

# ---- 1: signed ECDF difference maps, data vs model and null, same colour scale -------------------
rows = [("data_a.jsonl.1000.npz", "data vs model (seed 1000)"), ("null_0.jsonl.0.npz", "null: reference vs realisation 0")]
fig, axs = plt.subplots(2, 4, figsize=(13, 7.2), sharex=True, sharey=True, constrained_layout=True)
lim = 0.03
for r, (f, lab) in enumerate(rows):
    z = np.load(f"{H}/{f}")
    ref = z["ref"]
    for c, k in enumerate(("bb", "bf", "rb", "rf")):
        ax = axs[r, c]
        d = z["d_" + k]
        o = np.argsort(np.abs(d))                     # largest |d| drawn last
        sc = ax.scatter(ref[o, 0], ref[o, 1], c=d[o], s=0.6, cmap=DIV, vmin=-lim, vmax=lim, rasterized=True,
                        linewidths=0)
        if r == 0:
            ax.set_title(NAMES[k], color=INK, fontsize=10)
        ax.text(0.97, 0.04, f"95th pct |diff| {np.percentile(np.abs(d), 95):.4f}", transform=ax.transAxes,
                ha="right", va="bottom", color=INK, fontsize=8.5)
        ax.grid(color=GRID, lw=0.5)
        ax.set_xlim(0.35, 0.95)
        if c == 0:
            ax.set_ylabel(f"{lab}\nF606W")
        if r == 1:
            ax.set_xlabel("F606W - F814W")
axs[0, 0].invert_yaxis()
cb = fig.colorbar(sc, ax=axs, shrink=0.8, pad=0.01)
cb.set_label("ECDF(reference) - ECDF(sCMD) at each reference star")
fig.savefig(f"{H}/orientation_maps.png", dpi=150)
plt.close(fig)

# ---- 2: data mean vs its own null (mean, 3 sigma) per orientation, 95th percentile -------------
T = pd.read_csv(f"{H}/summary_table.csv")
sel = [("p95", "bb"), ("p95", "bf"), ("p95", "rb"), ("p95", "rf"), ("p95", "max4(stat)"), ("p95", "pointmax4"),
       ("p99", "Ying(bf,renorm,mask0)")]
lab = ["bluer & brighter\n(pipeline)", "bluer & fainter\n(Ying orientation)", "redder & brighter",
       "redder & fainter", "max of the 4\nstatistics (FF)", "per-star max\nof 4, then p95",
       "Ying full recipe\n(renorm, mask 0, p99)"]
fig, ax = plt.subplots(figsize=(10, 4.2), constrained_layout=True)
for j, (st, orient) in enumerate(sel):
    row = T[(T.stat == st) & (T.orientation == orient)].iloc[0]
    ax.plot([j, j], [row.null_mean - 3 * row.null_std, row.null_3sig], color="#eb6834", lw=2, solid_capstyle="round")
    ax.plot(j, row.null_mean, "o", color="#eb6834", ms=8, mec="#fcfcfb", mew=2)
    ax.plot(j, row.data_mean, "o", color="#2a78d6", ms=8, mec="#fcfcfb", mew=2)
    ax.errorbar(j, row.data_mean, yerr=row.data_std, color="#2a78d6", lw=2, capsize=0)
    ax.text(j + 0.12, row.data_mean, f"{row.sep:.1f} sigma", color=INK, va="center", fontsize=8)
ax.plot([], [], "o", color="#eb6834", label="null (split ASTs, 40 realisations): mean, line spans +-3 sigma")
ax.plot([], [], "o", color="#2a78d6", label="data vs model (12 realisations): mean +- std")
ax.set_xticks(range(len(sel)))
ax.set_xticklabels(lab, fontsize=8)
ax.set_ylabel("2D KS statistic")
ax.set_ylim(0, None)
ax.grid(axis="y", color=GRID, lw=0.5)
ax.legend(frameon=False, loc="upper left", fontsize=8)
ax.set_title("Best cell (run 10099, 11.4 Gyr, mu 13.38, E(V-I) 0.015): no orientation brings the fit inside its own null",
             color=INK, loc="left")
fig.savefig(f"{H}/orientation_summary.png", dpi=150)
plt.close(fig)
print("written", f"{H}/orientation_maps.png", f"{H}/orientation_summary.png")
