"""Golden-value regression: EEP tracks and isochrones must not drift.

The other self-checks test structure and known analytic answers; this one
pins the actual numbers.  Synthetic inputs (generated on the fly, so no data
files are needed) go through ``make_eep`` and ``make_iso``, and a sample of the
results is compared with ``tests/data/golden.json``.  Those reference values
were produced by the version of the code validated against the Fortran
binaries (see the README), so any change to EEP placement, the distance
metric, phases or the isochrone interpolation shows up here.

If you change the numerics *on purpose*, regenerate the reference with::

    PYEEP_REGEN_GOLDEN=1 python3 tests/test_golden.py

and say why in the commit message.
"""

from __future__ import annotations

import io
import json
import os
import sys
import tempfile
from contextlib import redirect_stdout

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

import make_eep  # noqa: E402
import test_make_eep  # noqa: E402
import test_pipeline  # noqa: E402
from pyeep.iso import read_eep  # noqa: E402

GOLDEN = os.path.join(HERE, 'data', 'golden.json')
REGEN = os.environ.get('PYEEP_REGEN_GOLDEN') == '1'
RTOL = 1e-12     # results are bit-reproducible; this only absorbs libm noise
EVERY = 25       # keep every 25th row so the reference file stays small


def _sample(rows: np.ndarray) -> list:
    return rows[::EVERY].tolist() + [rows[-1].tolist()]


def _eep_case():
    """The phased synthetic track of test_make_eep (exercises phases too)."""
    with tempfile.TemporaryDirectory() as root, redirect_stdout(io.StringIO()):
        test_make_eep._build_case(root)
        cwd = os.getcwd()
        try:
            os.chdir(root)
            assert make_eep.run('input.example') == 0
        finally:
            os.chdir(cwd)
        t = read_eep(os.path.join(root, 'eeps', '00100M.track.eep'))
    rows = np.vstack([t.tr, t.phase[None, :]]).T
    return {'eeps': t.eep.tolist(), 'rows': _sample(rows)}


def _pipeline_case(iso_nml):
    with redirect_stdout(io.StringIO()):
        _, tracks, isos = test_pipeline._pipeline(iso_nml + '\n   log_age_delta = 1.0')
    return {
        'tracks': {f"{t.initial_mass:.2f}": {'eeps': t.eep.tolist(),
                                             'rows': _sample(t.tr.T)}
                   for t in tracks[::3]},
        'isochrones': [{'age': iso.age, 'eeps': iso.eep.tolist(),
                        'rows': _sample(iso.data)} for iso in isos],
    }


CASES = {
    'eep_phased_1Msun': _eep_case,
    'pipeline_linear': lambda: _pipeline_case('do_linear_interpolation = .true.'),
    'pipeline_steffen': lambda: _pipeline_case('do_smooth = .false.\n   do_PAV = .false.'),
    'pipeline_defaults': lambda: _pipeline_case(''),
}


def _compare(path, got, want):
    if isinstance(want, dict):
        assert set(got) == set(want), f"{path}: keys {sorted(got)} != {sorted(want)}"
        for k in want:
            _compare(f"{path}.{k}", got[k], want[k])
    elif isinstance(want, list) and want and isinstance(want[0], dict):
        assert len(got) == len(want), f"{path}: {len(got)} items != {len(want)}"
        for i, (g, w) in enumerate(zip(got, want)):
            _compare(f"{path}[{i}]", g, w)
    elif isinstance(want, list) and want and isinstance(want[0], int):
        assert got == want, f"{path}: {got} != {want}"
    elif isinstance(want, list):
        g, w = np.array(got, dtype=float), np.array(want, dtype=float)
        assert g.shape == w.shape, f"{path}: shape {g.shape} != {w.shape}"
        bad = ~np.isclose(g, w, rtol=RTOL, atol=0.0, equal_nan=True)
        if bad.any():
            i = tuple(np.argwhere(bad)[0])
            raise AssertionError(f"{path}{list(i)}: {g[i]!r} != golden {w[i]!r} "
                                 f"({int(bad.sum())} values differ)")
    else:
        assert got == want, f"{path}: {got!r} != {want!r}"


def _check(name):
    got = CASES[name]()
    if REGEN:
        data = {}
        if os.path.exists(GOLDEN):
            with open(GOLDEN) as fh:
                data = json.load(fh)
        data[name] = got
        os.makedirs(os.path.dirname(GOLDEN), exist_ok=True)
        with open(GOLDEN, 'w') as fh:
            json.dump(data, fh, indent=0, sort_keys=True)
        print(f"regenerated {name}")
        return
    with open(GOLDEN) as fh:
        want = json.load(fh)[name]
    _compare(name, got, want)
    print(f"{name}: matches golden values")


def test_golden_eep_track():
    _check('eep_phased_1Msun')


def test_golden_pipeline_linear():
    _check('pipeline_linear')


def test_golden_pipeline_steffen():
    _check('pipeline_steffen')


def test_golden_pipeline_defaults():
    _check('pipeline_defaults')


if __name__ == '__main__':
    for n in CASES:
        _check(n)
    print('OK')
