#!/usr/bin/env python3
"""DSEP3 vs DSEP4 (both Monte Carlo ON) 47 Tuc isochrones, same var sets.

    compare_iso_d3_d4.py FIRST_RUN LAST_RUN [--file pyeep|example] [--tag NAME] [--d4 DIR]

Reads scmd_tb/tuc_iso/<run> (DSEP3, dsepX_fixed) and scmd_tb/tuc_iso_d4/<run>
(DSEP4), both built with the same pyeep recipe (8-16 Gyr in 0.1 Gyr steps); --d4 DIR
reads the DSEP4 isochrones from DIR instead (e.g. scmd_tb/tuc_iso_d4v2, stopped at the RGB tip).
--file pyeep reads the raw pyeep output (isochrones_pyeep.txt); --file example
(default) reads it after patch_lower_ms.py (isochrones_example.txt). Turnoff,
RGB and age-shift numbers are left blank for an isochrone that never reaches
the giant branch (max log L < 1). Writes into this directory:

  iso_overlay_<tag>.png   one panel per run: turnoff / subgiant / lower RGB at
                          10, 12 and 14 Gyr, DSEP3 solid, DSEP4 dotted
  iso_diffs_<tag>.png     DSEP4 - DSEP3 at 12 Gyr per run: main-sequence Teff at
                          log L = -1, turnoff Teff and L, RGB Teff at log L = 1.5,
                          and the age shift that matches the turnoff luminosity
  iso_diffs_<tag>.csv     the same numbers
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D3 = "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso/{run}/isos/isochrones_{file}.txt"
D4 = "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso_d4/{run}/isos/isochrones_{file}.txt"
COLS = ["EEP", "isochrone_age_yr", "initial_mass", "star_mass", "star_mdot", "he_core_mass",
        "c_core_mass", "log_L", "log_LH", "log_LHe", "log_Teff", "log_R", "log_g", "surface_h1",
        "surface_he3", "surface_he4", "surface_c12", "surface_o16", "log_center_T",
        "log_center_Rho", "center_gamma", "center_h1", "center_he4", "center_c12"]

# reference palette (dataviz skill), validated light-mode pair
C3, C4 = "#2a78d6", "#eb6834"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 9,
})


def read_iso(path):
    rows = [l.split() for l in open(path).read().split("\n")[11:]
            if l.strip() and not l.lstrip().startswith("#")]
    return pd.DataFrame(np.array(rows, dtype=float), columns=COLS)


def block(df, age):
    b = df[np.isclose(df.isochrone_age_yr, age, rtol=0, atol=1e6)]
    return b.sort_values("initial_mass")


def reaches_rgb(b):
    return len(b) > 2 and b.log_L.max() > 1.0


def turnoff(b):
    """Bluest point before the subgiant branch (log g > 3.5)."""
    pre = b[b.log_g > 3.5]
    k = pre.log_Teff.idxmax()
    return pre.loc[k, "log_Teff"], pre.loc[k, "log_L"], pre.loc[k, "initial_mass"]


def ms_teff(b, logl=-1.0):
    """log Teff on the unevolved main sequence at a given log L."""
    ms = b[(b.center_h1 > 0.3) & (b.log_L < 0.0)].sort_values("log_L")
    return np.interp(logl, ms.log_L, ms.log_Teff) if len(ms) > 2 and ms.log_L.min() < logl else np.nan


def rgb_teff(b, logl=1.5):
    """log Teff on the red giant branch at a given log L."""
    rgb = b[(b.center_h1 < 1e-6) & (b.log_g < 3.5) & (b.log_L > 0.8)].sort_values("log_L")
    return np.interp(logl, rgb.log_L, rgb.log_Teff) if len(rgb) > 2 else np.nan


def age_shift(d3, d4, age=12e9):
    """Age (Gyr) at which DSEP3's turnoff L equals DSEP4's turnoff L at `age`, minus `age`."""
    ages = np.sort(d3.isochrone_age_yr.unique())
    tol3 = np.array([turnoff(block(d3, a))[1] for a in ages])
    target = turnoff(block(d4, age))[1]
    o = np.argsort(tol3)  # turnoff L falls with age; interp needs increasing x
    return (np.interp(target, tol3[o], ages[o]) - age) / 1e9


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("first", type=int)
    ap.add_argument("last", type=int)
    ap.add_argument("--file", choices=("pyeep", "example"), default="example")
    ap.add_argument("--tag")
    ap.add_argument("--d4", help="directory of the DSEP4 runs (default tuc_iso_d4)")
    a = ap.parse_args()
    tag = a.tag or f"{a.file}_{a.first}_{a.last}"
    d3path = lambda r: D3.format(run=r, file=a.file)
    d4tmpl = D4 if a.d4 is None else os.path.join(a.d4, "{run}/isos/isochrones_{file}.txt")
    d4path = lambda r: d4tmpl.format(run=r, file=a.file)
    runs = [r for r in range(a.first, a.last + 1)
            if os.path.isfile(d3path(r)) and os.path.isfile(d4path(r))]
    skipped = sorted(set(range(a.first, a.last + 1)) - set(runs))
    if skipped:
        print("skipped (missing DSEP3 or DSEP4 isochrone):", skipped)

    ages = (10e9, 12e9, 14e9)
    ncol = 4
    nrow = int(np.ceil(len(runs) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(3.3 * ncol, 3.0 * nrow), squeeze=False)
    rec = []
    for i, r in enumerate(runs):
        d3, d4 = read_iso(d3path(r)), read_iso(d4path(r))
        ax = axs[i // ncol][i % ncol]
        for age in ages:
            b3, b4 = block(d3, age), block(d4, age)
            ax.plot(b3.log_Teff, b3.log_L, "-", color=C3, lw=1.4)
            ax.plot(b4.log_Teff, b4.log_L, ":", color=C4, lw=1.8)
            if not reaches_rgb(b3):
                continue
            t, l, _ = turnoff(b3)
            ax.annotate(f"{age / 1e9:.0f}", (t, l), xytext=(5, -2), textcoords="offset points",
                        color=INK2, fontsize=7, va="top")
        ax.set_xlim(3.84, 3.62)
        ax.set_ylim(-0.6, 1.8)
        n4 = sum(reaches_rgb(block(d4, age)) for age in ages)
        note = "" if n4 == len(ages) else f"   DSEP4 reaches RGB at {n4}/{len(ages)} ages"
        ax.set_title(f"run {r}{note}", loc="left", fontsize=8 if note else 9, color=INK)
        if i % ncol == 0:
            ax.set_ylabel("log L/L$_\\odot$")
        if i // ncol == nrow - 1 or i + ncol >= len(runs):
            ax.set_xlabel("log T$_{eff}$")
        b3, b4 = block(d3, 12e9), block(d4, 12e9)
        row = dict(run=r, d_ms_logTeff_at_logL_m1=ms_teff(b4) - ms_teff(b3),
                   d4_reaches_rgb_12gyr=reaches_rgb(b4))
        if reaches_rgb(b3) and reaches_rgb(b4):
            t3, l3, m3 = turnoff(b3)
            t4, l4, m4 = turnoff(b4)
            row.update(d_TO_logTeff=t4 - t3, d_TO_logL=l4 - l3, TO_mass_d3=m3, TO_mass_d4=m4,
                       d_RGB_logTeff_at_logL_1p5=rgb_teff(b4) - rgb_teff(b3),
                       age_shift_Gyr=age_shift(d3, d4))
        rec.append(row)
    for j in range(len(runs), nrow * ncol):
        axs[j // ncol][j % ncol].axis("off")
    h3, = axs[0][0].plot([], [], "-", color=C3, lw=1.4, label="DSEP3 (dsepX_fixed), Monte Carlo ON")
    h4, = axs[0][0].plot([], [], ":", color=C4, lw=1.8, label="DSEP4, Monte Carlo ON")
    fig.legend(handles=[h3, h4], loc="upper right", frameon=False, ncol=2, fontsize=9)
    fig.suptitle(f"47 Tuc isochrones ({a.file}), same var sets (runs {a.first}-{a.last}): 10, 12, 14 Gyr",
                 x=0.01, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(os.path.join(HERE, f"iso_overlay_{tag}.png"), dpi=160)
    plt.close(fig)

    df = pd.DataFrame(rec)
    for c in ("d_TO_logTeff", "d_TO_logL", "TO_mass_d3", "TO_mass_d4",
              "d_RGB_logTeff_at_logL_1p5", "age_shift_Gyr"):
        if c not in df:
            df[c] = np.nan
    df["age_shift_pct"] = 100 * df.age_shift_Gyr / 12.0
    df.to_csv(os.path.join(HERE, f"iso_diffs_{tag}.csv"), index=False, float_format="%.6f")

    panels = [("d_ms_logTeff_at_logL_m1", "MS log T$_{eff}$ at log L = -1"),
              ("d_TO_logTeff", "turnoff log T$_{eff}$"),
              ("d_TO_logL", "turnoff log L"),
              ("d_RGB_logTeff_at_logL_1p5", "RGB log T$_{eff}$ at log L = 1.5"),
              ("age_shift_pct", "age shift matching turnoff L (%)")]
    fig, axs = plt.subplots(1, len(panels), figsize=(2.6 * len(panels), 0.32 * len(df) + 1.6),
                            sharey=True)
    y = np.arange(len(df))
    for ax, (c, label) in zip(axs, panels):
        ax.axvline(0, color=INK2, lw=0.8)
        ax.plot(df[c], y, "o", color=C4, ms=5, mec=SURF, mew=1)
        med, n = df[c].median(), int(df[c].notna().sum())
        if n:
            ax.axvline(med, color=C4, lw=1, ls=":")
        ax.set_title(f"{label}\nmedian {med:+.4f} (n={n})", fontsize=8.5, loc="left")
        ax.grid(axis="y", visible=False)
    axs[0].set_yticks(y)
    axs[0].set_yticklabels([str(r) for r in df.run])
    axs[0].invert_yaxis()
    fig.suptitle(f"DSEP4 - DSEP3 at 12 Gyr, per run ({a.file} isochrones; blank = DSEP4 isochrone has no turnoff/RGB)",
                 x=0.01, ha="left", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(os.path.join(HERE, f"iso_diffs_{tag}.png"), dpi=160)
    plt.close(fig)
    pd.set_option("display.width", 200)
    print(df.round(5).to_string(index=False))


if __name__ == "__main__":
    main()
