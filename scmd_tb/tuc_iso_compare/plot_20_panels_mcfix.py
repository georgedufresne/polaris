#!/usr/bin/env python3
"""d3_d4_20_panels_hrd.png / _cmd.png redrawn with DSEP4 isochrones that read their Monte Carlo inputs.

The original figure (this folder, 2026-11-02) used scmd_tb/tuc_iso_d4v2, whose DSEP4 runs never read their var files:
~/nml_dsep4_mc set LMonteCarlo only for [alpha/Fe] = 0, and every 47 Tuc DSEP4 var set has 0.2 or 0.4 (no
'triple-alpha luminosity reached 1%' in their logs). Here DSEP4 = scmd_tb/tuc_iso_d4v2_mcfix: the same runs rerun with
LMonteCarlo on (tuc_dsep4_batch/run_tuc_dsep4_mcfix.batch) and rebuilt with the same v2 isochrone builder. DSEP3 is
unchanged (scmd_tb/tuc_iso, dsepX_fixed, which reads its var files). Panels, axes, colours and layout are those of
plot_d3_d4_isochrones.twenty_panels (imported); only the DSEP4 path, the title and the output folder (mcfix/) differ.
    plot_20_panels_mcfix.py [--with-old]   (--with-old also draws the old physics-off DSEP4 in grey)"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt     # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plot_d3_d4_isochrones as P   # noqa: E402

D4FIX = "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso_d4v2_mcfix/{run}/isos/isochrones_example.txt"
RUNS = [10000, 10002] + list(range(10004, 10022))      # the 20 panels of the original figure
OUT = os.path.join(HERE, "mcfix")
OLD = "#9a9a96"


def panels(runs, data, cmd, with_old):
    fig, axs = plt.subplots(4, 5, figsize=(15, 13), sharex=True, sharey=True)
    for ax, r in zip(axs.flat, runs):
        if with_old:
            for age in P.AGES:
                b = P.block(data["DSEP4 old"][r], age)
                x, y = (b.col, b.V) if cmd else (b.log_Teff, b.log_L)
                ax.plot(x, y, color=OLD, ls=P.LS[age], lw=1.0, alpha=0.6, zorder=1)
        # semi-transparent, DSEP3 a little wider underneath, so both show where they overlap
        for name, lw in (("DSEP3", 2.0), ("DSEP4", 1.2)):
            for age in P.AGES:
                b = P.block(data[name][r], age)
                x, y = (b.col, b.V) if cmd else (b.log_Teff, b.log_L)
                ax.plot(x, y, color=P.C[name], ls=P.LS[age], lw=lw, alpha=0.65, zorder=2)
        if cmd:
            P.shade_window(ax)
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
    h = [plt.Line2D([], [], color=P.C["DSEP3"], lw=2, label="DSEP3"),
         plt.Line2D([], [], color=P.C["DSEP4"], lw=2, label="DSEP4, physics on")]
    if with_old:
        h.append(plt.Line2D([], [], color=OLD, lw=2, label="DSEP4 old (physics off)"))
    h += [plt.Line2D([], [], color=P.INK2, lw=1.4, ls=P.LS[a], label=f"{a / 1e9:.0f} Gyr") for a in P.AGES]
    if cmd:
        h.append(plt.Rectangle((0, 0), 1, 1, color=P.BAND, label=f"2D KS window at $\\mu$ = {P.MU}"))
    fig.legend(handles=h, loc="upper right", ncol=len(h), frameon=False, fontsize=9)
    what = "CMD (MIST colours)" if cmd else "HR diagram"
    fig.suptitle(f"47 Tuc DSEP3 vs DSEP4, both reading their Monte Carlo inputs, 8 and 12 Gyr, one panel per run: {what}",
                 x=0.01, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    os.makedirs(OUT, exist_ok=True)
    out = os.path.join(OUT, f"d3_d4_20_panels_{'cmd' if cmd else 'hrd'}{'_with_old' if with_old else ''}.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


if __name__ == "__main__":
    runs = [r for r in RUNS if os.path.isfile(P.D3.format(run=r)) and os.path.isfile(D4FIX.format(run=r))]
    print(f"{len(runs)} of {len(RUNS)} runs have both DSEP3 and the fixed DSEP4 isochrones: {runs}")
    data = {"DSEP3": {r: P.read(P.D3.format(run=r)) for r in runs}, "DSEP4": {r: P.read(D4FIX.format(run=r)) for r in runs},
            "DSEP4 old": {r: P.read(P.D4.format(run=r)) for r in runs}}
    for with_old in (False, True):
        for cmd in (False, True):
            print(panels(runs, data, cmd, with_old))
