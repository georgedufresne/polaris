#!/usr/bin/env python3
"""NGC 6791 isochrone grid from W v1.4 (CONF1DENCE-W fcfac3c, frozen copy in this folder): the isochrone_W.py v1.4
recipe with the age table built once per run. A copy of ../emul_W/gen_W_iso.py (W v1.2, October), changed only where
the recipe changed.

isochrone_W.age_inversion() rebuilds the (EEP x mass) log-age table for every target age, but the table depends only on
the physics: the age enters as `tgt`. Here it is built once per run (same chunking) and, per age, the v1.4 root
selection (descending roots, EEP walked down from the tip with a mass cap, tol = 5: a root up to 5 mass steps above the
cap is kept), the 44 bisection halvings and the tip-anchored monotonic mass are copied from age_inversion, and the
position is isochrone_W.central_mode (the flow's mode). The EEP grid is isochrone_W.eep_grid (pyeep 1.4, EEP 202-1055 =
ZAMS to the RGB tip). `--validate RUN` checks that this gives the same isochrone as isochrone_W.w_isochrone.

Kept from the October run, so the pipeline is the same: ages 4.0-12.0 Gyr in 0.1, masses 0.10-2.5 Msun on 1200 points,
the var-file -> label mapping (21 values, C12+He4 = 1.0), Z from gs98getz, the MIST-format output.
Conductive opacity: W v1.4 stores it as a perturbation in [-1, 1] (default 0, the DSEP4 Blouin20 convention), while
these varfeh files hold a Hubbard-style multiplier (cocoef ~1 +/- 0.2). `--cond var` (default) passes the var-file value
as the October run did; `--cond zero` uses 0.0 (the mean Blouin opacity).
Inputs are copies on scratch: the repo snapshot, var/ and gs98getz in this folder.
--z-factor F (default 1, i.e. off): W's training grid computed Z with Martin's gs98getz, exactly 1.121 x ours, so W's
[Fe/H] label means Z = 1.121 x getz([Fe/H]). With F = 1.121 the flow and the age head are given
[Fe/H]_W = getz^-1(getz([Fe/H]) / F) at the run's [a/Fe], so W models the same Z as our DSEP3/OG runs of that var file
(emul_W14/feh_convention). The output header keeps the var file's [Fe/H] and Z = getz([Fe/H])."""
import argparse
import os
import subprocess
import sys
import time

import numpy as np

S = "/dartfs-hpc/scratch/f007058/emul_W14"
sys.path.insert(0, f"{S}/CONF1DENCE-W_fcfac3c")
from isochrone_W import central_mode, eep_grid, w_isochrone      # noqa: E402
from predict_W import WEmulator, AgeHead                         # noqa: E402

ISO_COLUMNS = ["EEP", "isochrone_age_yr", "initial_mass", "star_mass", "star_mdot",
               "he_core_mass", "c_core_mass", "log_L", "log_LH", "log_LHe",
               "log_Teff", "log_R", "log_g", "surface_h1", "surface_he3",
               "surface_he4", "surface_c12", "surface_o16", "log_center_T",
               "log_center_Rho", "center_gamma", "center_h1", "center_he4",
               "center_c12"]


def age_table(ah, p22, mlo, mhi, nm, chunk=256):
    """log age on the (EEP x mass) grid, as in isochrone_W.age_inversion (same grid and chunking)."""
    eg = eep_grid(ah)
    mg = np.linspace(mlo, mhi, nm)
    neep = len(eg)
    LA = np.concatenate([ah.log_age(np.repeat(eg[i:i + chunk], nm), np.tile(mg, len(eg[i:i + chunk])), p22).reshape(-1, nm)
                         for i in range(0, neep, chunk)])
    return eg, mg, LA


def invert(ah, p22, eg, mg, LA, age_gyr, mlo, mhi, tol=5, nbis=44):
    """Root selection + bisection + monotonic mass, copied from isochrone_W.age_inversion (v1.4)."""
    tgt = np.log10(age_gyr * 1e9)
    F = LA - tgt
    neep, nm = F.shape
    desc = (F[:, :-1] > 0) & (F[:, 1:] <= 0)
    sel = np.full(neep, -1, int); cap = nm - 2
    for r in range(neep - 1, -1, -1):
        c0 = np.nonzero(desc[r])[0]; c = c0[c0 <= cap]
        if len(c) == 0:
            a = c0[(c0 > cap) & (c0 <= cap + tol)]
            if len(a): sel[r] = a[0]
            continue
        sel[r] = c[-1]; cap = c[-1]
    rows = np.nonzero(sel >= 0)[0]
    if len(rows) == 0:
        return np.array([]), np.array([])
    j = sel[rows]; e = eg[rows]; lo = mg[j].copy(); hi = mg[j + 1].copy(); flo = F[rows, j]
    for _ in range(nbis):
        mid = 0.5 * (lo + hi); fm = ah.log_age(e, mid, p22) - tgt
        same = np.sign(fm) == np.sign(flo)
        lo = np.where(same, mid, lo); hi = np.where(same, hi, mid); flo = np.where(same, fm, flo)
    m = 0.5 * (lo + hi)
    m = np.maximum.accumulate(m) if tol == 0 else np.minimum.accumulate(m[::-1])[::-1]
    ok = (m > mlo + 1e-6) & (m < mhi - 1e-6)
    return m[ok], e[ok]


def phys_for(run, vardir, lh, default, cond):
    """22 physics in label_header order; C12+He4 (deboercoef) = 1.0 as upstream (as the October run)."""
    v = [float(l.split()[0]) for l in open(f'{vardir}/varfeh+031.{run}') if l.strip()]
    p = list(default)
    for nm, x in zip(lh[:22], [*v[0:16], 1.0, *v[16:21]]):
        p[lh.index(nm)] = x
    if cond == "zero":
        p[lh.index("Conductive Opacity")] = 0.0
    return p


def getz(feh, afe):
    return float(subprocess.run([f'{S}/gs98getz', f'{feh:.6f}', f'{afe:.6f}'], capture_output=True, text=True).stdout)


def feh_for_z(z, afe, lo=-2.0, hi=1.5):
    """Inverse of getz in [Fe/H] at fixed [a/Fe] (bisection, 50 halvings)."""
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if getz(mid, afe) < z else (lo, mid)
    return 0.5 * (lo + hi)


def write_mist(blocks, path, feh, y, z, afe):
    rows = []
    for age, eep, mass, o in blocks:
        n = len(eep)
        R = np.zeros((n, 24))
        R[:, 0] = eep
        R[:, 1] = age * 1e9
        R[:, 2] = mass
        R[:, 3] = mass
        R[:, 5] = o['he_core_mass']
        R[:, 7] = o['log_L']
        R[:, 10] = o['log_Teff']
        R[:, 11] = o['log_R']
        R[:, 12] = o['log_g']
        R[:, 22] = o['center_he4']
        rows.append(R)
    M = np.vstack(rows)
    npts = len(M) // max(len(blocks), 1)
    hdr = ("# MIST version number  = W14     \n"
           "# MESA revision number =    15140\n"
           "# " + "-" * 86 + "\n"
           "#  Yinit        Zinit   [Fe/H]   [a/Fe]  v/vcrit\n"
           f"# {y:.4f}  {z:.5E}     {feh:.2f}     {afe:.2f}     0.00\n"
           "# " + "-" * 86 + "\n"
           f"# number of isochrones = {len(blocks):5d}\n"
           "# " + "-" * 86 + "\n"
           f"# number of EEPs, cols = {npts:5d}   24\n"
           "# " + "  ".join(f"{i+1:>28d}" for i in range(24)) + "\n"
           "# " + "  ".join(f"{c:>28s}" for c in ISO_COLUMNS) + "\n")
    with open(path, 'w') as fh:
        fh.write(hdr)
        np.savetxt(fh, M, fmt='%28.16E')
    return len(M)


def validate(W, ah, lh, default, a):
    """Same isochrone as isochrone_W.w_isochrone? And how much does the conductive-opacity choice move it?"""
    run = a.validate
    res = {}
    for cond in ("var", "zero"):
        phys = phys_for(run, a.vardir, lh, default, cond)
        p22 = np.asarray(phys[:22], float)[:len(ah.theta)]
        t0 = time.time()
        eg, mg, LA = age_table(ah, p22, a.mass_low, a.mass_high, a.nmass)
        t1 = time.time()
        for age in (4.0, 8.0, 8.3, 12.0):
            m, e = invert(ah, p22, eg, mg, LA, age, a.mass_low, a.mass_high)
            o = central_mode(W, e, m, phys)
            res[(cond, age)] = (e, m, o)
            if cond == "var":
                ref = w_isochrone(W, ah, phys, age, mlo=a.mass_low, mhi=a.mass_high, nm=a.nmass)
                same = (len(ref["eep"]) == len(e) and np.array_equal(ref["eep"], e) and np.array_equal(ref["mass"], m)
                        and np.array_equal(ref["log_Teff"], o["log_Teff"]) and np.array_equal(ref["log_L"], o["log_L"]))
                print(f"run {run} {age:4.1f} Gyr: {len(e)} points, EEP {e.min():.1f}-{e.max():.1f}; identical to "
                      f"isochrone_W.w_isochrone: {same}", flush=True)
        print(f"  [cond={cond}] age table {t1 - t0:.1f} s; 4 ages {time.time() - t1:.1f} s", flush=True)
    cv = phys_for(run, a.vardir, lh, default, "var")[lh.index("Conductive Opacity")]
    for age in (4.0, 8.0, 8.3, 12.0):
        e0, m0, o0 = res[("var", age)]
        e1, m1, o1 = res[("zero", age)]
        dT = np.interp(e0, e1, o1["log_Teff"]) - o0["log_Teff"]
        dL = np.interp(e0, e1, o1["log_L"]) - o0["log_L"]
        k = o0["log_L"] < 1.5                                     # MS, turnoff, SGB, lower RGB
        print(f"conductive opacity 0.0 vs var value {cv:.3f}, {age:4.1f} Gyr: log L < 1.5: max |d log Teff| "
              f"{np.abs(dT[k]).max():.5f}, max |d log L| {np.abs(dL[k]).max():.5f}; whole isochrone: "
              f"{np.abs(dT).max():.5f}, {np.abs(dL).max():.5f}", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('runs', nargs='*')
    ap.add_argument('--vardir', default=f'{S}/var')
    ap.add_argument('--outdir', default='/dartfs-hpc/scratch/f007058/Isochrones_W14')
    ap.add_argument('--nmass', type=int, default=1200)
    ap.add_argument('--mass-low', type=float, default=0.10)
    ap.add_argument('--mass-high', type=float, default=2.5)
    ap.add_argument('--cond', choices=["var", "zero"], default="var")
    ap.add_argument('--device', default='cpu')
    ap.add_argument('--validate', type=int, default=None, help="run number: compare with isochrone_W.w_isochrone")
    ap.add_argument('--z-factor', type=float, default=1.0, help="W [Fe/H] for Z = getz/F (1.121: our Z convention)")
    a = ap.parse_args()

    W = WEmulator(device=a.device)
    ah = AgeHead(device=a.device)
    lh = W.fmt['label_header']
    default = W.fmt['label_default']
    if a.validate is not None:
        return validate(W, ah, lh, default, a)

    ages = np.round(np.arange(4.0, 12.0001, 0.1), 4)
    print(f"W v1.4 recipe (hoisted table) on {a.device}: {len(ages)} ages, EEP grid {len(eep_grid(ah))}, "
          f"nmass={a.nmass}, cond={a.cond}, z-factor={a.z_factor}", flush=True)
    for run in a.runs:
        vf = f'{a.vardir}/varfeh+031.{run}'
        if not os.path.exists(vf):
            print(f"  {run}: no var file", flush=True)
            continue
        out = f'{a.outdir}/{run}/isos/isochrones_example.txt'
        if os.path.exists(out) and os.path.getsize(out) > 1e6:
            print(f"  {run}: exists", flush=True)
            continue
        t0 = time.time()
        phys = phys_for(run, a.vardir, lh, default, a.cond)
        p22 = np.asarray(phys[:22], float)[:len(ah.theta)]
        feh, yy, afe = phys[0], phys[1], phys[2]
        z = float(subprocess.run([f'{S}/gs98getz', str(feh), str(afe)], capture_output=True, text=True).stdout.strip() or 0)
        if a.z_factor != 1.0:
            fw = feh_for_z(getz(feh, afe) / a.z_factor, afe)
            phys[lh.index("FeH")] = fw
            p22 = np.asarray(phys[:22], float)[:len(ah.theta)]
            print(f"  {run}: [Fe/H] {feh:+.6f} -> W [Fe/H] {fw:+.6f} (Z {z:.6e}; W grid Z {getz(fw, afe) * a.z_factor:.6e})", flush=True)
        eg, mg, LA = age_table(ah, p22, a.mass_low, a.mass_high, a.nmass)
        blocks = []
        for age in ages:
            m, e = invert(ah, p22, eg, mg, LA, float(age), a.mass_low, a.mass_high)
            if len(e) < 10:
                continue
            o = central_mode(W, e, m, phys)
            blocks.append((float(age), e, m, o))
        if not blocks:
            print(f"  {run}: no isochrones", flush=True)
            continue
        os.makedirs(os.path.dirname(out), exist_ok=True)
        tmp = out + ".tmp"
        n = write_mist(blocks, tmp, feh, yy, z, afe)
        os.replace(tmp, out)
        print(f"  {run} -> {len(blocks)} ages, {n:,} rows, {os.path.getsize(out)/1e6:.1f} MB, {time.time() - t0:.0f} s", flush=True)


if __name__ == '__main__':
    main()
