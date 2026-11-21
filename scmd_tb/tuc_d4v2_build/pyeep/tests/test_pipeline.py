"""End-to-end: synthetic MESA history files -> make_eep -> make_iso.

Builds a small grid of synthetic tracks whose ages scale exactly as
``age = age_1Msun * M**-2.5``.  Shifting log(age) by a constant does not
change the EEP distance metric, so every track gets the same EEPs, and at a
given EEP the age--mass relation is known analytically.  The test therefore
checks both halves of the pipeline:

* EEPs: every track gets an ``.eep`` file with identical primary EEPs and the
  scaled ages at each EEP;
* isochrones: the initial masses make_iso recovers match the analytic
  inverse of the age--mass relation.

Run directly (``python tests/test_pipeline.py``) or via ``run_tests.py``.
"""

from __future__ import annotations

import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PY_ROOT = os.path.dirname(HERE)
sys.path.insert(0, PY_ROOT)
sys.path.insert(0, HERE)

import make_eep  # noqa: E402
import make_iso  # noqa: E402
from make_synthetic_track import COLUMNS, _synthetic_arrays, write_history_columns_list  # noqa: E402
from pyeep.iso import read_eep  # noqa: E402
from test_make_eep import INPUT_EEP, INPUT_FORMAT  # noqa: E402

MASSES = [0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40]
SLOPE = 2.5
AGES = [3.0e8, 6.0e8]   # inside the synthetic tracks (they end at 0.4-1.5 Gyr)

NML = """&eep_controls
   eep_input_file = 'input.eep'
   do_phases = .false.
   center_gamma_limit = 100
   center_carbon_limit = 0.05
   log_center_T_limit = 8.5
   high_mass_limit = 10.0
   very_low_mass_limit = 0.7
   weight_center_rho_T_by_Xc = .false.
   Tc_scale = 1d-2
   Rhoc_scale = 1d-2
   age_scale = 1d0
   Teff_scale = 2.0d0
   logL_scale = 0.125d0
/

&iso_controls
   {iso}
/
"""


def _write_history(path, m):
    """Same layout as make_synthetic_track, ages scaled by M**-SLOPE."""
    data = _synthetic_arrays(m)
    data['star_age'] = data['star_age'] * m ** -SLOPE
    n = len(data['star_age'])
    num = lambda k: ' '.join(f"{i:>28d}" for i in range(1, k + 1))  # noqa: E731
    lines = [num(3),
             ' '.join(f"{h:>28}" for h in ('version_number', 'initial_mass', 'initial_z')),
             ' '.join(f"{v:>28}" for v in ('12000', f"{m:.6f}", '0.014')),
             '', num(len(COLUMNS)), ' '.join(f"{c:>28}" for c in COLUMNS)]
    lines += [' '.join(f"{data[c][j]:>28.16e}" for c in COLUMNS) for j in range(n)]
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


def _build(root, iso_nml):
    hist = os.path.join(root, 'hist')
    os.makedirs(hist)
    names = []
    for m in MASSES:
        name = f"m{int(round(m * 100)):04d}"
        _write_history(os.path.join(hist, name), m)
        names.append(name)
    write_history_columns_list(os.path.join(root, 'cols.list'))
    for fname, text in (('input.eep', INPUT_EEP), ('input.format', INPUT_FORMAT),
                        ('input.nml', NML.format(iso=iso_nml))):
        with open(os.path.join(root, fname), 'w') as fh:
            fh.write(text)
    with open(os.path.join(root, 'input.example'), 'w') as fh:
        fh.write(f"#v\npipe\n#Y Z FeH aFe v\n 0.2700 1.40000e-02 0.00 0.00 0.0\n"
                 f"#dirs\n{hist}\n{root}/eeps\n{root}/isos\n#cols\n{root}/cols.list\n"
                 f"#tracks\n{len(names)}\n" + ''.join(n + '\n' for n in names)
                 + f"#isochrones\npipe.iso\nlist\nlinear\n{len(AGES)}\n"
                 + ''.join(f"{a:.6e}\n" for a in AGES))
    return names


def _pipeline(iso_nml):
    """Run make_eep then make_iso in a scratch directory; return what they made."""
    with tempfile.TemporaryDirectory() as root:
        names = _build(root, iso_nml)
        cwd = os.getcwd()
        try:
            os.chdir(root)
            assert make_eep.run('input.example') == 0
            isos = make_iso.make_isochrones('input.example')
        finally:
            os.chdir(cwd)
        tracks = [read_eep(os.path.join(root, 'eeps', n + '.eep')) for n in names]
        wrote_iso = os.path.exists(os.path.join(root, 'isos', 'pipe.iso'))
    return wrote_iso, tracks, isos


def test_pipeline_eeps():
    wrote_iso, tracks, _ = _pipeline('do_linear_interpolation = .true.')
    ref = tracks[0]
    assert ref.neep >= 4 and ref.eep[0] == 1
    i_age = ref.col_names.index('star_age')
    for m, t in zip(MASSES, tracks):
        assert t.initial_mass == m
        assert t.eep.tolist() == ref.eep.tolist(), "same primary EEPs for every mass"
        # at every EEP the age is the 0.8 Msun age scaled by (m/0.8)**-SLOPE
        scale = (m / MASSES[0]) ** -SLOPE
        assert np.allclose(t.tr[i_age], ref.tr[i_age] * scale, rtol=1e-9)
    assert wrote_iso
    print(f"pipeline EEPs: {len(tracks)} tracks, primary EEPs {ref.eep.tolist()}")


def test_pipeline_isochrones():
    worst = {}
    for label, iso_nml, tol in (
            ('linear', 'do_linear_interpolation = .true.', 5e-3),
            ('Steffen', 'do_smooth = .false.\n   do_PAV = .false.', 5e-4)):
        _, tracks, isos = _pipeline(iso_nml + '\n   log_age_delta = 1.0')
        ref = tracks[0]
        age_ref = ref.tr[ref.col_names.index('star_age')]       # 0.8 Msun ages
        for iso, age in zip(isos, AGES):
            assert iso.neep > 20, iso.neep
            mi = iso.column('initial_mass')
            assert np.all(np.diff(mi) >= -1e-12), "initial mass must rise along the isochrone"
            assert mi.min() >= MASSES[0] - 1e-12 and mi.max() <= MASSES[-1] + 1e-12
            # analytic: age = age_ref(eep) * (m/0.8)**-SLOPE
            m_true = MASSES[0] * (age / age_ref[iso.eep - 1]) ** (-1 / SLOPE)
            err = np.max(np.abs(mi - m_true))
            assert err < tol, (label, age, err)
            worst[label] = max(worst.get(label, 0.0), err)
            assert np.allclose(iso.column('star_mass'), mi, atol=tol)
    print("pipeline isochrones: worst initial-mass error "
          + ", ".join(f"{k} {v:.1e} Msun" for k, v in worst.items()))


if __name__ == '__main__':
    test_pipeline_eeps()
    test_pipeline_isochrones()
    print("OK")
