"""Self-checks for the make_iso port (no MESA data or Fortran needed).

* the vectorised kernels (Steffen interpolation, binary search, PAV, output
  formatting) against straightforward scalar ports of the Fortran/MESA code,
  bit-for-bit, on random inputs;
* an end-to-end run on a synthetic grid of ``.eep`` tracks, written with
  pyeep's own ``write_track``, checking the recovered initial masses against
  the analytic age--mass relation the grid was built from.

Run directly (``python tests/test_make_iso.py``) or under pytest.
"""

from __future__ import annotations

import math
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
PY_ROOT = os.path.dirname(HERE)
sys.path.insert(0, PY_ROOT)

import make_iso  # noqa: E402
from pyeep.iso import (  # noqa: E402
    IsoControls, IsochroneGrid, MakeIsoError, _binary_search, format_rows,
    pav, read_eep, steffen_eval,
)
from pyeep.support import (  # noqa: E402
    Column, EEPState, Track, STAR_LOW_MASS, _fmt_es32_16e3, binary_search,
    write_track,
)


# ---------------------------------------------------------------------------
# scalar references, transcribed directly from the Fortran
# ---------------------------------------------------------------------------

def ref_pm_value(x, f, xv):
    """mk_pmcub + do_interp_values, full-array, one value."""
    n = len(x)
    c2, c3, c4 = [0.0] * n, [0.0] * n, [0.0] * n
    if n == 2:
        c2[0] = (f[1] - f[0]) / (x[1] - x[0])
    else:
        h = [x[i + 1] - x[i] for i in range(n - 1)]
        s = [(f[i + 1] - f[i]) / h[i] for i in range(n - 1)]
        for i in range(1, n - 1):
            p = (s[i - 1] * h[i] + s[i] * h[i - 1]) / (h[i - 1] + h[i])
            c2[i] = ((math.copysign(1.0, s[i - 1]) + math.copysign(1.0, s[i]))
                     * min(abs(s[i - 1]), abs(s[i]), 0.5 * abs(p)))
        p = s[0] * (1 + h[0] / (h[0] + h[1])) - s[1] * h[0] / (h[0] + h[1])
        c2[0] = 0.0 if p * s[0] <= 0 else (2 * s[0] if abs(p) > 2 * abs(s[0]) else p)
        p = (s[n - 2] * (1 + h[n - 2] / (h[n - 2] + h[n - 3]))
             - s[n - 3] * h[n - 2] / (h[n - 2] + h[n - 3]))
        c2[n - 1] = 0.0 if p * s[n - 2] <= 0 else (
            2 * s[n - 2] if abs(p) > 2 * abs(s[n - 2]) else p)
        for i in range(n - 1):
            c3[i] = (3 * s[i] - 2 * c2[i] - c2[i + 1]) / h[i]
            c4[i] = (c2[i] + c2[i + 1] - 2 * s[i]) / (h[i] * h[i])
    inc = x[0] < x[n - 1]
    if inc:
        xv = x[n - 1] if xv > x[n - 1] else (x[0] if xv < x[0] else xv)
    else:
        xv = x[n - 1] if xv < x[n - 1] else (x[0] if xv > x[0] else xv)
    k, xk, xk1 = 0, x[0], x[1]
    while (inc and xv > xk1) or ((not inc) and xv < xk1):
        k += 1
        if k >= n - 1:
            k -= 1
            xv = xk1
            break
        xk, xk1 = xk1, x[k + 1]
    d = xv - xk
    return f[k] + d * (c2[k] + d * (c3[k] + d * c4[k]))


def test_steffen_matches_reference_bitwise():
    rng = np.random.default_rng(42)
    cases = 0
    for n in (2, 3, 4, 5, 8, 20):
        for _ in range(300):
            x = np.cumsum(rng.uniform(0.01, 1.0, n))
            if rng.random() < 0.5:
                x = x[::-1].copy()                    # decreasing abscissae
            f = rng.normal(size=n)
            if rng.random() < 0.3:
                f = np.sort(f)
            lo, hi = min(x[0], x[-1]), max(x[0], x[-1])
            xv = rng.uniform(lo - 0.1, hi + 0.1)
            if rng.random() < 0.1:
                xv = float(x[rng.integers(n)])        # exactly on a node
            got = steffen_eval(x[None, :], np.array([n]), np.array([xv]),
                               lambda p: f[p])[0]
            want = ref_pm_value(x.tolist(), f.tolist(), xv)
            assert got == want or (math.isnan(got) and math.isnan(want)), \
                (n, x, f, xv, got, want)
            cases += 1
    print(f"steffen_eval: {cases} random cases bit-identical")


def test_binary_search_matches_scalar():
    rng = np.random.default_rng(7)
    for _ in range(3000):
        n = int(rng.integers(1, 12))
        v = rng.normal(size=n)
        if rng.random() < 0.7:
            v = np.sort(v)[::(1 if rng.random() < 0.5 else -1)].copy()
        val = float(rng.choice(v)) if rng.random() < 0.2 else float(rng.normal())
        got = int(_binary_search(v[None, :], np.array([n]), val)[0])
        assert got == binary_search(n, v, 1, val), (v, val)
    print("binary_search: 3000 random cases identical")


def test_pav():
    y = pav(np.array([0.70, 0.72, 0.74, 0.76, 0.78, 0.49]))
    assert np.all(np.diff(y) >= 0)
    assert np.allclose(y, 0.6983333333333334)
    assert pav(np.array([1.0, 2.0, 3.0])).tolist() == [1.0, 2.0, 3.0]
    y = pav(np.array([3.0, 1.0, 2.0, 0.0, 5.0]))
    assert np.all(np.diff(y) >= 0) and y[-1] == 5.0
    print("pav: OK")


def test_format_rows():
    vals = np.array([[1.0, -2.5e-7, 0.0, -0.0, 1.0e100, 3.3e-120, 4.0e9]])
    fast = format_rows(np.array([7]), vals[:, [0, 1, 2, 3, 6]])
    slow = f"{7:5d}" + ''.join(_fmt_es32_16e3(v) for v in vals[0, [0, 1, 2, 3, 6]]) + '\n'
    assert fast == slow, (fast, slow)
    wide = format_rows(np.array([7]), vals)
    assert wide == f"{7:5d}" + ''.join(_fmt_es32_16e3(v) for v in vals[0]) + '\n'
    assert '1.0000000000000000E+100' in wide and '3.2999999999999997E-120' in wide
    assert all(len(r) == 5 + 32 * vals.shape[1] for r in wide.splitlines())
    nan = format_rows(np.array([1]), np.array([[np.nan, np.inf]]))
    assert nan == '    1' + 'NaN'.rjust(32) + 'Infinity'.rjust(32) + '\n'
    print("format_rows: OK")


# ---------------------------------------------------------------------------
# end to end on a synthetic grid
# ---------------------------------------------------------------------------

COLS = ['star_age', 'star_mass', 'log_L', 'log_Teff', 'center_h1']
NEEP = 120
MASSES = [0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.20, 1.30, 1.40, 1.50]


def lifetime(m):
    return 1.0e10 * m ** -2.5


def eep_age(m, j):
    """star_age at EEP j (1-based): fraction of the lifetime, monotone in j."""
    return lifetime(m) * (0.02 + 0.98 * (j - 1) / (NEEP - 1)) ** 1.5


def write_grid(root, kink=None):
    eep_dir = os.path.join(root, 'eeps')
    os.makedirs(eep_dir)
    state = EEPState()
    names = []
    for m in MASSES:
        j = np.arange(1, NEEP + 1)
        age = eep_age(m, j)
        if kink is not None:          # stretch one EEP's age in every track
            age[kink - 1] *= 1.3
        frac = (j - 1) / (NEEP - 1)
        t = Track()
        t.version_string = 'synth'
        t.MESA_revision_number = 15140
        t.initial_mass = m
        t.initial_Y, t.initial_Z, t.Fe_div_H = 0.27, 0.014, 0.0
        t.cols = [Column(c) for c in COLS]
        t.ncol = len(COLS)
        t.ntrack = NEEP
        t.eep = np.array([1, 40, 80, NEEP])
        t.neep = t.eep.size
        t.star_type = STAR_LOW_MASS
        t.tr = np.vstack([age, np.full(NEEP, m),
                          4.5 * np.log10(m) + 2.0 * frac,
                          3.76 + 0.6 * np.log10(m) - 0.15 * frac,
                          0.7 * (1 - frac)])
        name = f"m{int(round(m * 100)):04d}"
        t.filename = os.path.join(eep_dir, name + '.eep')
        write_track(t, state)
        names.append(name)
    with open(os.path.join(root, 'cols.list'), 'w') as fh:
        fh.write('\n'.join(COLS) + '\n')
    return eep_dir, names


def write_input(root, eep_dir, names, ages, style=''):
    p = os.path.join(root, 'input.example')
    with open(p, 'w') as fh:
        fh.write("#v\nsynth\n#Y Z\n 0.2700 1.40000e-02 0.00 0.00 0.0\n#dirs\n"
                 f"{root}/hist\n{eep_dir}\n{root}/isos\n#cols\n"
                 f"{root}/cols.list\n#tracks\n{len(names)}\n"
                 + ''.join(n + '\n' for n in names)
                 + "#isochrones\nout.iso\nlist\nlinear\n"
                 f"{len(ages)}\n" + ''.join(f"{a:.6e}\n" for a in ages)
                 + (style + '\n' if style else ''))
    return p


def test_read_eep_round_trip():
    with tempfile.TemporaryDirectory() as root:
        eep_dir, names = write_grid(root)
        t = read_eep(os.path.join(eep_dir, names[4] + '.eep'))
        assert t.initial_mass == 1.0 and t.ntrack == NEEP and t.ncol == len(COLS)
        assert t.eep.tolist() == [1, 40, 80, NEEP]
        assert t.star_type == STAR_LOW_MASS and not t.has_phase
        assert t.version_string == 'synth   ' and t.MESA_revision_number == 15140
        assert t.tr[0, -1] == eep_age(1.0, NEEP)
    print("read_eep: OK")


def _run(root, nml, ages, style='', kink=None):
    eep_dir, names = write_grid(root, kink)
    inp = write_input(root, eep_dir, names, ages, style)
    with open(os.path.join(root, 'input.nml'), 'w') as fh:
        fh.write('&iso_controls\n' + nml + '/\n')
    cwd = os.getcwd()
    try:
        os.chdir(root)
        isos = make_iso.make_isochrones(inp)
    finally:
        os.chdir(cwd)
    return isos, os.path.join(root, 'isos', 'out.iso')


def test_end_to_end_recovers_masses():
    ages = [2.0e9, 4.0e9, 6.0e9]
    for nml, tol in (('do_linear_interpolation=.true.\n', 1e-2),
                     ('do_smooth=.false.\ndo_PAV=.false.\n', 5e-4),
                     ('', 1e-2)):
        with tempfile.TemporaryDirectory() as root:
            isos, out = _run(root, nml + 'log_age_delta=1.0\n', ages)
            assert os.path.exists(out)
            for iso, age in zip(isos, ages):
                assert iso.neep > 10, (nml, iso.neep)
                mi = iso.column('initial_mass')
                assert np.all(np.diff(mi) >= -1e-12), "masses must be monotone"
                # invert eep_age analytically: lifetime(m) = age / frac**1.5
                fr = (0.02 + 0.98 * (iso.eep - 1) / (NEEP - 1)) ** 1.5
                m_true = (age / fr / 1.0e10) ** (-1 / 2.5)
                assert np.max(np.abs(mi - m_true)) < tol, (nml, np.max(np.abs(mi - m_true)))
                # star_mass is initial mass in this grid, so it must agree
                assert np.allclose(iso.column('star_mass'), mi, atol=tol)
            with open(out) as fh:
                text = fh.read()
            assert text.count('# number of EEPs') == len(ages)
            assert '# number of isochrones =     3' in text
    print("end-to-end (linear / Steffen / defaults): masses recovered")


def test_pav_enforces_monotonic_mass():
    """A kink in one EEP's ages gives a mass reversal only PAV removes."""
    kink, age = 60, [4.0e9]
    base = 'do_linear_interpolation=.true.\ndo_smooth=.false.\nlog_age_delta=1.0\n'
    with tempfile.TemporaryDirectory() as root:
        raw, _ = _run(root, base + 'do_PAV=.false.\n', age, kink=kink)
    with tempfile.TemporaryDirectory() as root:
        fixed, _ = _run(root, base + 'do_PAV=.true.\n', age, kink=kink)
    m_raw = raw[0].column('initial_mass')
    m_fix = fixed[0].column('initial_mass')
    assert (np.diff(m_raw) < 0).any(), "kinked grid should reverse without PAV"
    assert np.all(np.diff(m_fix) >= 0), "PAV must leave initial mass monotone"
    assert raw[0].eep.tolist() == fixed[0].eep.tolist()
    changed = np.flatnonzero(m_raw != m_fix)
    assert changed.size >= 2 and kink in raw[0].eep[changed]
    # pooled points share one mean value, and the total is conserved
    assert np.allclose(m_fix[changed], m_fix[changed][0])
    assert np.isclose(m_fix[changed].sum(), m_raw[changed].sum())
    print(f"PAV: pooled EEPs {raw[0].eep[changed].tolist()}")


def test_double_and_errors():
    with tempfile.TemporaryDirectory() as root:
        isos, _ = _run(root, 'do_smooth=.false.\n', [4.0e9], 'double')
        assert isos[0].neep > 0
    with tempfile.TemporaryDirectory() as root:
        eep_dir, names = write_grid(root)
        names = [names[1], names[2], names[0]] + names[3:]   # 3rd < 1st
        inp = write_input(root, eep_dir, names, [4.0e9])
        with open(os.path.join(root, 'input.nml'), 'w') as fh:
            fh.write('&iso_controls\n/\n')
        cwd = os.getcwd()
        try:
            os.chdir(root)
            try:
                make_iso.make_isochrones(inp)
            except MakeIsoError as exc:
                assert 'out of order' in str(exc)
            else:
                raise AssertionError('expected masses-out-of-order error')
            os.remove('input.nml')
            try:
                make_iso.make_isochrones(inp)
            except MakeIsoError as exc:
                assert 'namelist' in str(exc)
            else:
                raise AssertionError('missing input.nml must be an error')
        finally:
            os.chdir(cwd)
    print("double-EEP mode and error paths: OK")


def test_parallel_batch_identical():
    with tempfile.TemporaryDirectory() as root:
        eep_dir, names = write_grid(root)
        with open(os.path.join(root, 'input.nml'), 'w') as fh:
            fh.write('&iso_controls\n/\n')
        inputs = []
        for i in range(4):
            sub = os.path.join(root, f'g{i}')
            os.makedirs(sub)
            p = write_input(sub, eep_dir, names, [2e9 + i * 1e9, 5e9])
            with open(p) as fh:
                txt = fh.read().replace(f'{sub}/cols.list', f'{root}/cols.list')
            with open(p, 'w') as fh:
                fh.write(txt)
            inputs.append(p)
        cwd = os.getcwd()
        try:
            os.chdir(root)
            assert make_iso.run(inputs, jobs=1) == 0
            serial = [open(os.path.join(root, f'g{i}', 'isos', 'out.iso')).read()
                      for i in range(4)]
            assert make_iso.run(inputs, jobs=3) == 0
            par = [open(os.path.join(root, f'g{i}', 'isos', 'out.iso')).read()
                   for i in range(4)]
        finally:
            os.chdir(cwd)
        assert serial == par
    print("parallel batch: identical to serial")


if __name__ == '__main__':
    test_steffen_matches_reference_bitwise()
    test_binary_search_matches_scalar()
    test_pav()
    test_format_rows()
    test_read_eep_round_trip()
    test_end_to_end_recovers_masses()
    test_pav_enforces_monotonic_mass()
    test_double_and_errors()
    test_parallel_batch_identical()
    print("OK")
