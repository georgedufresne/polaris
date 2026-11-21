"""Generate a small synthetic MESA-style history file for testing.

Real MESA history files are not shipped with the repository (and the Fortran
program requires a full MESA install to build), so this module fabricates a
physically-plausible-enough low-mass track that exercises the EEP finder:
central hydrogen depletes (PreMS -> ZAMS -> TAMS), then central helium burns
and central gamma rises so the star is classified low-mass and reaches the
white-dwarf cooling sequence.

The output file follows the MESA history layout expected by
``pyeep.support.read_history_file``:

    line 1: header column numbers
    line 2: header names          (includes version_number, initial_mass)
    line 3: header values
    line 4: (blank)
    line 5: data column numbers
    line 6: data column names
    line 7+: data rows
"""

from __future__ import annotations

import os
from typing import Dict, List

import numpy as np

COLUMNS = [
    'star_age', 'star_mass', 'log_LH', 'log_LHe', 'log_Teff', 'log_L',
    'log_g', 'log_center_T', 'log_center_Rho', 'center_h1', 'center_he4',
    'center_c12', 'center_gamma', 'surface_h1', 'he_core_mass', 'c_core_mass',
]


def _synthetic_arrays(initial_mass: float = 1.0, n: int = 400) -> Dict[str, np.ndarray]:
    idx = np.arange(n)

    # three phases: pre-MS/MS (H burn), RGB+He burn, WD cooling
    n_ms = n // 2
    n_he = n // 4

    star_age = np.geomspace(1.0e4, 1.2e10, n)
    star_mass = np.full(n, initial_mass)

    # central hydrogen: flat pre-MS, then depletes to ~0 by end of MS
    center_h1 = np.empty(n)
    n_pre = 10
    center_h1[:n_pre] = 0.70
    center_h1[n_pre:n_ms] = np.linspace(0.70, 1.0e-20, n_ms - n_pre)
    center_h1[n_ms:] = 1.0e-20

    # central helium: flat during MS, burns during He phase, then ~0
    center_he4 = np.empty(n)
    center_he4[:n_ms] = 0.28
    center_he4[n_ms:n_ms + n_he] = np.linspace(0.28, 1.0e-8, n_he)
    center_he4[n_ms + n_he:] = 1.0e-8

    # central carbon: builds up as helium burns
    center_c12 = np.empty(n)
    center_c12[:n_ms] = 1.0e-10
    center_c12[n_ms:n_ms + n_he] = np.linspace(1.0e-10, 0.5, n_he)
    center_c12[n_ms + n_he:] = 0.5

    # log g: rises to a maximum near ZAMS (~point n_pre) then declines
    log_g = np.empty(n)
    log_g[:n_pre] = np.linspace(3.5, 4.4, n_pre)
    log_g[n_pre:] = np.linspace(4.4, 1.0, n - n_pre)

    # log L: modest on MS, peaks at the RGB tip (end of H-shell / He flash)
    log_L = np.empty(n)
    log_L[:n_ms] = np.linspace(0.0, 0.6, n_ms)
    log_L[n_ms:n_ms + n_he] = np.concatenate([
        np.linspace(0.6, 3.4, n_he // 2),          # climb to RGB tip
        np.linspace(3.4, 1.6, n_he - n_he // 2),   # settle onto HB
    ])
    log_L[n_ms + n_he:] = np.linspace(1.6, -3.0, n - n_ms - n_he)

    # log Teff: cool at the RGB tip, hot on the WD cooling track
    log_Teff = np.empty(n)
    log_Teff[:n_ms] = np.linspace(3.76, 3.70, n_ms)
    log_Teff[n_ms:n_ms + n_he] = np.concatenate([
        np.linspace(3.70, 3.55, n_he // 2),        # coolest at RGB tip
        np.linspace(3.55, 3.85, n_he - n_he // 2),
    ])
    log_Teff[n_ms + n_he:] = np.linspace(3.85, 4.9, n - n_ms - n_he)

    # central temperature: rises overall, dips a little onto ZAHB
    log_center_T = np.empty(n)
    log_center_T[:n_ms] = np.linspace(7.1, 7.6, n_ms)
    log_center_T[n_ms:n_ms + n_he] = np.concatenate([
        np.linspace(7.6, 8.2, n_he // 2),
        np.linspace(8.2, 8.05, n_he - n_he // 2),  # slight dip = ZAHB
    ])
    log_center_T[n_ms + n_he:] = np.linspace(8.05, 7.5, n - n_ms - n_he)

    log_center_Rho = np.linspace(1.8, 7.0, n)

    # He/He-shell burning power
    log_LH = np.empty(n)
    log_LH[:n_ms] = np.linspace(-0.2, 0.6, n_ms)
    log_LH[n_ms:] = np.linspace(0.6, -6.0, n - n_ms)

    log_LHe = np.full(n, -30.0)
    log_LHe[n_ms:n_ms + n_he] = np.linspace(-2.0, 2.5, n_he)  # He ignition
    log_LHe[n_ms + n_he:] = np.linspace(2.5, -8.0, n - n_ms - n_he)

    # central electron degeneracy: crosses the WD threshold (100) at the end
    center_gamma = np.empty(n)
    center_gamma[:n_ms + n_he] = np.linspace(0.01, 5.0, n_ms + n_he)
    center_gamma[n_ms + n_he:] = np.linspace(5.0, 250.0, n - n_ms - n_he)

    surface_h1 = np.full(n, 0.70)

    # core masses grow through the giant phases
    he_core_mass = np.empty(n)
    he_core_mass[:n_ms] = np.linspace(0.0, 0.1, n_ms)
    he_core_mass[n_ms:] = np.linspace(0.1, 0.55, n - n_ms)
    c_core_mass = np.empty(n)
    c_core_mass[:n_ms + n_he] = np.linspace(0.0, 0.2, n_ms + n_he)
    c_core_mass[n_ms + n_he:] = np.linspace(0.2, 0.54, n - n_ms - n_he)

    return {
        'star_age': star_age,
        'star_mass': star_mass,
        'log_LH': log_LH,
        'log_LHe': log_LHe,
        'log_Teff': log_Teff,
        'log_L': log_L,
        'log_g': log_g,
        'log_center_T': log_center_T,
        'log_center_Rho': log_center_Rho,
        'center_h1': center_h1,
        'center_he4': center_he4,
        'center_c12': center_c12,
        'center_gamma': center_gamma,
        'surface_h1': surface_h1,
        'he_core_mass': he_core_mass,
        'c_core_mass': c_core_mass,
    }


def write_history_file(path: str, initial_mass: float = 1.0,
                       version_number: int = 12000, n: int = 400) -> None:
    data = _synthetic_arrays(initial_mass, n)
    ncol = len(COLUMNS)

    header_names = ['version_number', 'initial_mass', 'initial_z']
    header_vals = [str(version_number), f"{initial_mass:.6f}", "0.014"]

    def numrow(ncols: int) -> str:
        return ' '.join(f"{k:>28d}" for k in range(1, ncols + 1))

    lines: List[str] = []
    lines.append(numrow(len(header_names)))
    lines.append(' '.join(f"{h:>28}" for h in header_names))
    lines.append(' '.join(f"{v:>28}" for v in header_vals))
    lines.append('')
    lines.append(numrow(ncol))
    lines.append(' '.join(f"{c:>28}" for c in COLUMNS))
    for j in range(n):
        row = ' '.join(f"{data[c][j]:>28.16e}" for c in COLUMNS)
        lines.append(row)

    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


def write_history_columns_list(path: str) -> None:
    """Write a minimal history_columns.list containing the required columns."""
    with open(path, 'w') as fh:
        fh.write("! synthetic history columns for testing\n\n")
        for c in COLUMNS:
            fh.write(f"      {c}\n")


if __name__ == '__main__':
    import sys
    out = sys.argv[1] if len(sys.argv) > 1 else 'synthetic.data'
    write_history_file(out)
    print(f"wrote {out}")
