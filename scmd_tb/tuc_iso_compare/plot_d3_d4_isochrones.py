#!/usr/bin/env python3
"""DSEP3 vs DSEP4 47 Tuc isochrones at 8 and 12 Gyr, same Monte Carlo var sets.

DSEP3: scmd_tb/tuc_iso/<run> (dsepX_fixed, ALPHAC = 0); DSEP4: scmd_tb/tuc_iso_d4v2/<run>
(Monte Carlo ON, pyeep stopped at the RGB tip). Both isochrones_example.txt, the files the
sCMD grids read; CMD magnitudes from the grids' own MIST ACS/WFC table ([Fe/H] -0.75).
The shaded band is the 2D KS window (F606W 15.02-19.02) at mu = 13.39. Writes here:

  d3_d4_one_run_<run>.png   one run: HR diagram and CMD, 8 and 12 Gyr on the same axes
  d3_d4_20_runs.png         the first 20 runs with both, 8 and 12 Gyr, both models, same axes
  d3_d4_all_runs.png       every run with both: HR diagram and CMD, 8 and 12 Gyr
  d3_d4_differences.png     DSEP4 - DSEP3 colour at fixed M_F606W along the MS and RGB, every
                            run (median and 16-84% band), and the turnoff / subgiant shifts

    plot_d3_d4_isochrones.py [RUN]     (RUN for the one-run figure; default 10091)
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
B = "/dartfs-hpc/scratch/f007058/scmd_tb"
sys.path.insert(0, f"{B}/tuc_grid_v4/code")
from isochronetoolbox.scmd.model import BCTable  # noqa: E402

D3 = B + "/tuc_iso/{run}/isos/isochrones_example.txt"
D4 = B + "/tuc_iso_d4v2/{run}/isos/isochrones_example.txt"
COLS = ["EEP", "isochrone_age_yr", "initial_mass", "star_mass", "star_mdot", "he_core_mass",
        "c_core_mass", "log_L", "log_LH", "log_LHe", "log_Teff", "log_R", "log_g", "surface_h1",
        "surface_he3", "surface_he4", "surface_c12", "surface_o16", "log_center_T",
        "log_center_Rho", "center_gamma", "center_h1", "center_he4", "center_c12"]
AGES = (8e9, 12e9)
BC = BCTable(f"{B}/tuc_grid_v3/code/data/47tuc/bc_mist_acs_wfc_fehm075.npz")
MU, WIN = 13.39, (15.02, 19.02)

# validated light-mode pair (dataviz reference palette), recessive ink and grid
C = {"DSEP3": "#2a78d6", "DSEP4": "#eb6834"}
LS = {8e9: (0, (4, 2)), 12e9: "-"}
INK, INK2, MUTED, GRID, SURF, BAND = "#0b0b0b", "#52514e", "#8a8985", "#e4e3df", "#fcfcfb", "#efeee9"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 9,
    "legend.frameon": False, "legend.fontsize": 8,
})


def order(b):
    """One age block in evolutionary order: up to the turnoff by initial mass (the rows
    patch_lower_ms.py added have EEP = 0 and lie on the lower MS, where mass rises along the
    isochrone), after it by EEP (near the RGB tip the masses differ by ~1e-3 and are not
    monotonic, so sorting by mass scrambles the giant branch)."""
    real = b[b.EEP > 0]
    e_to = real.loc[real[real.log_g > 3.5].log_Teff.idxmax(), "EEP"]
    pre = b[(b.EEP == 0) | (b.EEP <= e_to)].sort_values("initial_mass")
    post = b[b.EEP > e_to].sort_values("EEP")
    return pd.concat([pre, post])


def read(path):
    """The 8 and 12 Gyr blocks of an isochrone file in evolutionary order, with absolute F606W, F814W."""
    d = pd.read_csv(path, sep=r"\s+", comment="#", header=None, names=COLS, engine="c")
    d = d[np.isclose(d.isochrone_age_yr, AGES[0], atol=1e6) | np.isclose(d.isochrone_age_yr, AGES[1], atol=1e6)].copy()
    d["age"] = np.where(np.isclose(d.isochrone_age_yr, AGES[0], atol=1e6), AGES[0], AGES[1])
    m = BC.absolute(d)
    d["V"], d["col"] = m["F606W"], m["F606W"] - m["F814W"]
    return pd.concat([order(d[d.age == a]) for a in AGES]).reset_index(drop=True)


def block(d, age):
    return d[d.age == age]


def turnoff(b):
    """Index of the bluest point before the subgiant branch (log g > 3.5)."""
    pre = b[b.log_g > 3.5]
    return pre.log_Teff.idxmax()


def branches(b):
    """(MS, RGB) parts of one age block, each sorted by M_F606W: MS below the turnoff; RGB
    above the subgiant branch (log g < 3.6, brighter than the turnoff by 1 mag) up to its
    brightest F606W -- beyond that the coolest giants turn fainter in F606W, so M_F606W
    stops being a coordinate along the branch (DSEP4 reaches the tip, DSEP3 stops near
    log L 2.6, so only DSEP4 shows it)."""
    k = turnoff(b)
    ms = b.loc[:k]
    ms = ms[ms.V > b.loc[k, "V"] + 0.3].sort_values("V")
    post = b.loc[k:]
    rgb = post[(post.log_g < 3.6) & (post.V < b.loc[k, "V"] - 1.0)]
    if len(rgb):
        rgb = rgb.iloc[:int(np.argmin(rgb.V.values)) + 1]
    return ms, rgb.sort_values("V")


def sgb_v(b):
    """M_F606W of the subgiant point 0.05 mag redder than the turnoff (as MakeTestCMD)."""
    k = turnoff(b)
    post = b.loc[k:]
    j = np.argmax(post.col.values >= b.loc[k, "col"] + 0.05)
    return post.V.values[j]


def shade_window(ax, cmd=True):
    if cmd:
        ax.axhspan(WIN[0] - MU, WIN[1] - MU, color=BAND, zorder=0, lw=0)


def label_ages(ax, b3, cmd):
    for age in AGES:
        b = block(b3, age)
        k = turnoff(b)
        x, y = (b.loc[k, "col"], b.loc[k, "V"]) if cmd else (b.loc[k, "log_Teff"], b.loc[k, "log_L"])
        ax.annotate(f"{age / 1e9:.0f} Gyr", (x, y), xytext=(-6, 0) if cmd else (6, 0), textcoords="offset points",
                    ha="right" if cmd else "left", va="center", fontsize=8, color=INK2)


def hrd_axes(ax):
    ax.set_xlim(3.84, 3.58)
    ax.set_ylim(-1.3, 2.7)
    ax.set_xlabel("log T$_{eff}$")
    ax.set_ylabel("log L/L$_\\odot$")


def cmd_axes(ax):
    ax.set_xlim(0.40, 1.15)
    ax.set_ylim(7.2, -1.6)
    ax.set_xlabel("F606W $-$ F814W")
    ax.set_ylabel("M$_{F606W}$")


def legend(ax, extra=()):
    h = [plt.Line2D([], [], color=C[m], lw=2, ls="-", label=m) for m in C]
    h += [plt.Line2D([], [], color=INK2, lw=1.4, ls=LS[a], label=f"{a / 1e9:.0f} Gyr") for a in AGES]
    h += list(extra)
    ax.legend(handles=h, loc="upper left")


def runs_with_both():
    return [r for r in range(10000, 11000) if os.path.isfile(D3.format(run=r)) and os.path.isfile(D4.format(run=r))]


def one_run(run):
    d3, d4 = read(D3.format(run=run)), read(D4.format(run=run))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 5.2))
    for name, d in (("DSEP3", d3), ("DSEP4", d4)):
        for age in AGES:
            b = block(d, age)
            a1.plot(b.log_Teff, b.log_L, color=C[name], ls=LS[age], lw=1.8)
            a2.plot(b.col, b.V, color=C[name], ls=LS[age], lw=1.8)
    shade_window(a2)
    hrd_axes(a1); cmd_axes(a2)
    label_ages(a1, d3, False); label_ages(a2, d3, True)
    legend(a1)
    legend(a2, [plt.Rectangle((0, 0), 1, 1, color=BAND, label=f"2D KS window at $\\mu$ = {MU}")])
    a1.set_title("HR diagram", loc="left"); a2.set_title("CMD (MIST ACS/WFC colours, [Fe/H] $-$0.75)", loc="left")
    v = dict((l.split()[1], float(l.split()[0])) for l in open(f"/dartfs-hpc/scratch/f007058/47Tuc_var/varfeh-077.{run}") if l.strip())
    fig.suptitle(f"47 Tuc run {run}: DSEP3 vs DSEP4 isochrones at 8 and 12 Gyr   "
                 f"([Fe/H] {v['FeH']:+.2f}, Y {v['Yprim']:.3f}, [$\\alpha$/Fe] {v['alphafe']:.1f}, "
                 f"$\\alpha_{{MLT}}$ {v['CMIXLA']:.2f}, KTTAU {int(v['KTTAU'])})", x=0.01, ha="left", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out = os.path.join(HERE, f"d3_d4_one_run_{run}.png")
    fig.savefig(out, dpi=170); plt.close(fig)
    return out


def twenty_runs(runs, data):
    """The first 20 runs with both models: 8 and 12 Gyr, DSEP3 and DSEP4, all on the same axes."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 6))
    for name in ("DSEP3", "DSEP4"):
        for r in runs:
            for age in AGES:
                b = block(data[name][r], age)
                a1.plot(b.log_Teff, b.log_L, color=C[name], ls=LS[age], lw=0.9, alpha=0.55)
                a2.plot(b.col, b.V, color=C[name], ls=LS[age], lw=0.9, alpha=0.55)
    hrd_axes(a1); cmd_axes(a2); shade_window(a2)
    label_ages(a1, data["DSEP3"][runs[0]], False); label_ages(a2, data["DSEP3"][runs[0]], True)
    legend(a1)
    legend(a2, [plt.Rectangle((0, 0), 1, 1, color=BAND, label=f"2D KS window at $\\mu$ = {MU}")])
    a1.set_title("HR diagram", loc="left"); a2.set_title("CMD (MIST ACS/WFC colours, [Fe/H] $-$0.75)", loc="left")
    fig.suptitle(f"47 Tuc: DSEP3 vs DSEP4 isochrones at 8 and 12 Gyr, {len(runs)} Monte Carlo runs "
                 f"({runs[0]}-{runs[-1]}, same var sets)", x=0.01, ha="left", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(HERE, "d3_d4_20_runs.png")
    fig.savefig(out, dpi=170); plt.close(fig)
    return out


def twenty_panels(runs, data, cmd):
    """One panel per run (first 20 with both models): 8 and 12 Gyr, DSEP3 and DSEP4."""
    fig, axs = plt.subplots(4, 5, figsize=(15, 13), sharex=True, sharey=True)
    for ax, r in zip(axs.flat, runs):
        for name in ("DSEP3", "DSEP4"):
            for age in AGES:
                b = block(data[name][r], age)
                x, y = (b.col, b.V) if cmd else (b.log_Teff, b.log_L)
                ax.plot(x, y, color=C[name], ls=LS[age], lw=1.3)
        if cmd:
            shade_window(ax)
        v = dict((l.split()[1], float(l.split()[0])) for l in open(f"/dartfs-hpc/scratch/f007058/47Tuc_var/varfeh-077.{r}") if l.strip())
        ax.set_title(f"run {r}  [Fe/H] {v['FeH']:+.2f}  Y {v['Yprim']:.3f}", loc="left", fontsize=8)
    for ax in axs.flat[len(runs):]:
        ax.axis("off")
    if cmd:
        axs[0, 0].set_xlim(0.40, 1.05); axs[0, 0].set_ylim(6.5, -0.5)
    else:
        axs[0, 0].set_xlim(3.82, 3.62); axs[0, 0].set_ylim(-0.8, 2.2)
    for ax in axs[-1, :]:
        ax.set_xlabel("F606W $-$ F814W" if cmd else "log T$_{eff}$")
    for ax in axs[:, 0]:
        ax.set_ylabel("M$_{F606W}$" if cmd else "log L/L$_\\odot$")
    h = [plt.Line2D([], [], color=C[m], lw=2, label=m) for m in C]
    h += [plt.Line2D([], [], color=INK2, lw=1.4, ls=LS[a], label=f"{a / 1e9:.0f} Gyr") for a in AGES]
    if cmd:
        h.append(plt.Rectangle((0, 0), 1, 1, color=BAND, label=f"2D KS window at $\\mu$ = {MU}"))
    fig.legend(handles=h, loc="upper right", ncol=len(h), frameon=False, fontsize=9)
    what = "CMD (MIST colours)" if cmd else "HR diagram"
    fig.suptitle(f"47 Tuc DSEP3 vs DSEP4, 8 and 12 Gyr, one panel per run: {what}",
                 x=0.01, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    out = os.path.join(HERE, f"d3_d4_20_panels_{'cmd' if cmd else 'hrd'}.png")
    fig.savefig(out, dpi=150); plt.close(fig)
    return out


def all_runs(runs, data):
    fig, axs = plt.subplots(2, 2, figsize=(10.5, 9.5))
    for j, age in enumerate(AGES):
        for name in ("DSEP3", "DSEP4"):
            for r in runs:
                b = block(data[name][r], age)
                axs[0, j].plot(b.log_Teff, b.log_L, color=C[name], ls="-", lw=0.5, alpha=0.22)
                axs[1, j].plot(b.col, b.V, color=C[name], ls="-", lw=0.5, alpha=0.22)
        hrd_axes(axs[0, j]); cmd_axes(axs[1, j]); shade_window(axs[1, j])
        axs[0, j].set_title(f"{age / 1e9:.0f} Gyr: HR diagram", loc="left")
        axs[1, j].set_title(f"{age / 1e9:.0f} Gyr: CMD", loc="left")
    h = [plt.Line2D([], [], color=C[m], lw=2, label=f"{m} ({len(runs)} runs)") for m in C]
    axs[0, 0].legend(handles=h, loc="upper left")
    axs[1, 0].legend(handles=[plt.Rectangle((0, 0), 1, 1, color=BAND, label=f"2D KS window at $\\mu$ = {MU}")], loc="upper left")
    fig.suptitle(f"47 Tuc: every Monte Carlo run with both a DSEP3 and a DSEP4 isochrone ({len(runs)} runs, same var sets)",
                 x=0.01, ha="left", fontsize=10.5)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    out = os.path.join(HERE, "d3_d4_all_runs.png")
    fig.savefig(out, dpi=170); plt.close(fig)
    return out


def differences(runs, data):
    grid_ms, grid_rgb = np.arange(2.2, 7.01, 0.05), np.arange(-1.0, 2.01, 0.05)
    rows, curves = [], {a: {"ms": [], "rgb": []} for a in AGES}
    for r in runs:
        for age in AGES:
            b3, b4 = block(data["DSEP3"][r], age), block(data["DSEP4"][r], age)
            (m3, g3), (m4, g4) = branches(b3), branches(b4)
            for part, (x3, x4, grid) in (("ms", (m3, m4, grid_ms)), ("rgb", (g3, g4, grid_rgb))):
                lo, hi = max(x3.V.min(), x4.V.min()), min(x3.V.max(), x4.V.max())
                y = np.interp(grid, x4.V, x4.col) - np.interp(grid, x3.V, x3.col)
                curves[age][part].append(np.where((grid >= lo) & (grid <= hi), y, np.nan))
            k3, k4 = turnoff(b3), turnoff(b4)
            rows.append(dict(run=r, age=age / 1e9, d_to_col=b4.loc[k4, "col"] - b3.loc[k3, "col"],
                             d_sgb_v=sgb_v(b4) - sgb_v(b3)))
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(HERE, "d3_d4_differences.csv"), index=False, float_format="%.5f")

    fig = plt.figure(figsize=(11, 8.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1])
    axm = fig.add_subplot(gs[0, :2]); axr = fig.add_subplot(gs[0, 2])
    for age in AGES:
        for ax, part, grid in ((axm, "ms", grid_ms), (axr, "rgb", grid_rgb)):
            Y = 1000 * np.array(curves[age][part])
            for y in Y:
                ax.plot(y, grid, color=INK2, lw=0.4, alpha=0.10 if age == AGES[1] else 0.0)
            med, p16, p84 = (np.nanpercentile(Y, q, axis=0) for q in (50, 16, 84))
            col = C["DSEP4"] if age == AGES[1] else C["DSEP3"]
            ax.fill_betweenx(grid, p16, p84, color=col, alpha=0.15, lw=0)
            ax.plot(med, grid, color=col, lw=2, ls=LS[age], label=f"{age / 1e9:.0f} Gyr: median, 16-84%")
    for ax, title in ((axm, "main sequence (below the turnoff)"), (axr, "red giant branch")):
        ax.axvline(0, color=INK2, lw=0.8)
        ax.axhspan(WIN[0] - MU, WIN[1] - MU, color=BAND, zorder=0, lw=0)
        ax.set_xlabel("DSEP4 $-$ DSEP3 colour at fixed M$_{F606W}$ [mmag]")
        ax.set_title(title, loc="left")
    axm.set_ylim(7.0, 2.2); axr.set_ylim(2.0, -1.0)
    axm.set_ylabel("M$_{F606W}$")
    axm.legend(loc="lower left")
    axm.annotate("thin lines: each run at 12 Gyr", (0.99, 0.02), xycoords="axes fraction", ha="right", fontsize=7.5, color=MUTED)

    for j, (q, lab) in enumerate((("d_to_col", "turnoff colour, DSEP4 $-$ DSEP3 [mmag]"),
                                   ("d_sgb_v", "subgiant M$_{F606W}$, DSEP4 $-$ DSEP3 [mmag]"))):
        ax = fig.add_subplot(gs[1, j])
        for age in AGES:
            v = 1000 * t[t.age == age / 1e9][q]
            col = C["DSEP4"] if age == AGES[1] else C["DSEP3"]
            ax.hist(v, bins=np.linspace(np.nanpercentile(1000 * t[q], 1) - 5, np.nanpercentile(1000 * t[q], 99) + 5, 30),
                    histtype="step", lw=2, color=col, ls=LS[age], label=f"{age / 1e9:.0f} Gyr: median {v.median():+.1f}")
        ax.axvline(0, color=INK2, lw=0.8)
        ax.set_xlabel(lab); ax.set_ylabel("runs")
        ax.legend(loc="upper left")
        ax.grid(axis="x", visible=False)
    ax = fig.add_subplot(gs[1, 2])
    for age in AGES:
        s = t[t.age == age / 1e9]
        col = C["DSEP4"] if age == AGES[1] else C["DSEP3"]
        ax.scatter(1000 * s.d_to_col, 1000 * s.d_sgb_v, s=14, color=col, edgecolor=SURF, lw=0.6,
                   marker="o" if age == AGES[1] else "s", label=f"{age / 1e9:.0f} Gyr")
    ax.axvline(0, color=INK2, lw=0.8); ax.axhline(0, color=INK2, lw=0.8)
    ax.set_xlabel("turnoff colour shift [mmag]"); ax.set_ylabel("subgiant M$_{F606W}$ shift [mmag]")
    ax.legend(loc="upper left")
    ax.set_title("per run", loc="left")
    fig.suptitle(f"47 Tuc: DSEP4 $-$ DSEP3 isochrones, {len(runs)} runs, 8 and 12 Gyr "
                 f"(negative colour = DSEP4 bluer; negative M = DSEP4 brighter; shaded = KS window at $\\mu$ = {MU})",
                 x=0.01, ha="left", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    out = os.path.join(HERE, "d3_d4_differences.png")
    fig.savefig(out, dpi=170); plt.close(fig)
    return out, t


if __name__ == "__main__":
    run = int(sys.argv[1]) if len(sys.argv) > 1 else 10091
    print(one_run(run))
    runs = runs_with_both()
    data = {"DSEP3": {r: read(D3.format(run=r)) for r in runs}, "DSEP4": {r: read(D4.format(run=r)) for r in runs}}
    print(twenty_runs(runs[:20], data))
    print(all_runs(runs, data))
    out, t = differences(runs, data)
    print(out)
    for (age, q), v in t.melt(id_vars=["run", "age"], value_vars=["d_to_col", "d_sgb_v"]).groupby(["age", "variable"]).value:
        p = 1000 * np.nanpercentile(v, [16, 50, 84])
        print(f"{age:4.0f} Gyr {q:9s} [mmag]: median {p[1]:+6.1f}  16-84% {p[0]:+6.1f} .. {p[2]:+6.1f}  (n {v.notna().sum()})")
