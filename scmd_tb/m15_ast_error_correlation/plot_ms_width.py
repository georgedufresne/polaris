#!/usr/bin/env python3
"""ms_width.png: M15 main-sequence colour width (core, 1.4826 x MAD about the ridge) against F606W, observed and as
predicted by the two artificial-star error models (ms_width.csv from ast_error_correlation.py). Palette: categorical
slots 1-3, validated (validate_palette.js, light: all checks pass; aqua below 3:1 contrast, so every line is
labelled directly)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt     # noqa: E402
import pandas as pd                 # noqa: E402

H = "/dartfs-hpc/scratch/f007058/scmd_tb/m15_ast_error_correlation"
BLUE, ORANGE, AQUA, INK, MUTED, GRID, BG = "#2a78d6", "#eb6834", "#1baf7a", "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": BG, "axes.facecolor": BG})
d = pd.read_csv(f"{H}/ms_width.csv")
x = 0.5 * (d.f606w_lo + d.f606w_hi)
series = [("independent", ORANGE, "independent draws per filter (MakeTestCMD)"),
          ("observed", BLUE, "observed main sequence"),
          ("joint", AQUA, "both errors from one artificial star (ours)")]
fig, ax = plt.subplots(figsize=(8.0, 4.9))
for col, c, name in series:
    ax.plot(x, d[col], color=c, lw=2, marker="o", ms=8, mec=BG, mew=2, label=name)
    ax.annotate(name, (x.iloc[-1], d[col].iloc[-1]), xytext=(10, 0), textcoords="offset points", va="center",
                color=INK, fontsize=8.5)
ax.set_xlim(x.min() - 0.08, x.max() + 0.62)
ax.set_ylim(0, 0.045)
ax.set_xticks([19.4, 19.6, 19.8, 20.0, 20.2, 20.4, 20.6])
ax.set_xlabel("F606W")
ax.set_ylabel("colour width of the main sequence (mag)")
ax.grid(axis="y", color=GRID, lw=0.8)
ax.legend(frameon=False, loc="upper left", fontsize=8.5)
ax.set_title("M15: independent per-filter errors predict a wider main sequence than observed", color=INK,
             fontsize=10.5, loc="left")
ax.text(0.0, -0.2, "Observed width$^2$ = intrinsic$^2$ + error$^2$, so a model above the observed line is excluded.\n"
        "Artificial-star errors in the two filters: Spearman +0.68, Pearson +0.84 (clipped).", transform=ax.transAxes,
        color=MUTED, fontsize=8, va="top")
fig.tight_layout()
fig.savefig(f"{H}/ms_width.png", dpi=150)
print(f"{H}/ms_width.png")
