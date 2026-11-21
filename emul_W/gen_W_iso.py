#!/usr/bin/env python3
"""NGC 6791 isochrone grid from W -- official isochrone_W.py recipe, age table hoisted.

isochrone_W.age_inversion() rebuilds the (neep x nmass) age table for every target
age, but that table depends only on the physics; the age enters solely as `tgt`.
Over 81 ages per run that is 81x redundant (~27M age-head calls). Here the table
is built ONCE per run and the per-age work is just the root selection plus
bisection, which is what the official routine does anyway. The root selection
(descending roots, walk EEP down from the tip with a mass cap), the 44 bisection
halvings and np.maximum.accumulate are copied verbatim, and the position comes
from isochrone_W.central_mode -- so the output is identical, only faster.
"""
import argparse
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, '/dartfs-hpc/rc/home/8/f007058/CONF1DENCE-W')
from isochrone_W import central_mode, EEP_LO, EEP_HI           # noqa: E402
from predict_W import WEmulator, AgeHead                        # noqa: E402

ISO_COLUMNS = ["EEP", "isochrone_age_yr", "initial_mass", "star_mass", "star_mdot",
               "he_core_mass", "c_core_mass", "log_L", "log_LH", "log_LHe",
               "log_Teff", "log_R", "log_g", "surface_h1", "surface_he3",
               "surface_he4", "surface_c12", "surface_o16", "log_center_T",
               "log_center_Rho", "center_gamma", "center_h1", "center_he4",
               "center_c12"]


def age_table(ah, p22, mlo, mhi, nm, neep):
    """log age on the (EEP x mass) grid -- independent of the target age."""
    eg = np.linspace(EEP_LO, EEP_HI, neep)
    mg = np.linspace(mlo, mhi, nm)
    E = np.repeat(eg, nm)
    M = np.tile(mg, neep)
    return eg, mg, ah.log_age(E, M, p22).reshape(neep, nm)


def invert(ah, p22, eg, mg, LA, age_gyr, mlo, mhi, nbis=44):
    """Root selection + bisection, verbatim from isochrone_W.age_inversion."""
    tgt = np.log10(age_gyr * 1e9)
    F = LA - tgt
    neep, nm = F.shape
    desc = (F[:, :-1] > 0) & (F[:, 1:] <= 0)
    sel = np.full(neep, -1, int)
    cap = nm - 2
    for r in range(neep - 1, -1, -1):
        c = np.nonzero(desc[r])[0]
        c = c[c <= cap]
        if len(c) == 0:
            continue
        sel[r] = c[-1]
        cap = c[-1]
    rows = np.nonzero(sel >= 0)[0]
    if len(rows) == 0:
        return np.array([]), np.array([])
    j = sel[rows]
    e = eg[rows]
    lo = mg[j].copy()
    hi = mg[j + 1].copy()
    flo = F[rows, j]
    for _ in range(nbis):
        mid = 0.5 * (lo + hi)
        fm = ah.log_age(e, mid, p22) - tgt
        same = np.sign(fm) == np.sign(flo)
        lo = np.where(same, mid, lo)
        hi = np.where(same, hi, mid)
        flo = np.where(same, fm, flo)
    m = np.maximum.accumulate(0.5 * (lo + hi))
    ok = (m > mlo + 1e-6) & (m < mhi - 1e-6)
    return m[ok], e[ok]


def phys_for(run, vardir, lh, default):
    """22 physics in label_header order; C12+He4 (deboercoef) = 1.0 as upstream."""
    v = [float(l.split()[0]) for l in open(f'{vardir}/varfeh+031.{run}') if l.strip()]
    p = list(default)
    for nm, x in zip(lh[:22], [*v[0:16], 1.0, *v[16:21]]):
        p[lh.index(nm)] = x
    return p


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
    hdr = ("# MIST version number  = W       \n"
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('runs', nargs='+')
    ap.add_argument('--vardir', default='/dartfs-hpc/rc/home/8/f007058/mcdsep/MC_var_noco')
    ap.add_argument('--outdir', default='/dartfs-hpc/scratch/f007058/Isochrones_W')
    ap.add_argument('--neep', type=int, default=280)
    ap.add_argument('--nmass', type=int, default=1200)
    ap.add_argument('--mass-low', type=float, default=0.10)
    ap.add_argument('--mass-high', type=float, default=2.5)
    ap.add_argument('--device', default='cpu')
    a = ap.parse_args()

    W = WEmulator(device=a.device)
    ah = AgeHead(device=a.device)
    lh = W.fmt['label_header']
    default = W.fmt['label_default']
    ages = np.round(np.arange(4.0, 12.0001, 0.1), 4)
    print(f"W official recipe (hoisted table) on {a.device}: {len(ages)} ages, "
          f"neep={a.neep}, nmass={a.nmass}", flush=True)

    for run in a.runs:
        vf = f'{a.vardir}/varfeh+031.{run}'
        if not os.path.exists(vf):
            print(f"  {run}: no var file", flush=True)
            continue
        out = f'{a.outdir}/{run}/isos/isochrones_example.txt'
        if os.path.exists(out) and os.path.getsize(out) > 1e6:
            print(f"  {run}: exists", flush=True)
            continue
        phys = phys_for(run, a.vardir, lh, default)
        p22 = np.asarray(phys[:22], float)[:len(ah.theta)]
        feh, yy, afe = phys[0], phys[1], phys[2]
        z = float(subprocess.run(
            ['/dartfs-hpc/rc/home/8/f007058/mcdsep/run/gs98getz', str(feh), str(afe)],
            capture_output=True, text=True).stdout.strip() or 0)
        eg, mg, LA = age_table(ah, p22, a.mass_low, a.mass_high, a.nmass, a.neep)
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
        n = write_mist(blocks, out, feh, yy, z, afe)
        print(f"  {run} -> {len(blocks)} ages, {n:,} rows, "
              f"{os.path.getsize(out)/1e6:.1f} MB", flush=True)


if __name__ == '__main__':
    main()
