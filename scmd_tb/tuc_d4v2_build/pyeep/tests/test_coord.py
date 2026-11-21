"""Per-model EEP coordinate (pyeep.coord.eep_coordinate) on the synthetic track.

Run directly (``python tests/test_coord.py``) or under pytest.
"""

from __future__ import annotations

import copy
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.path.insert(0, HERE)

from pyeep.coord import eep_coordinate  # noqa: E402
from pyeep.eep import secondary_eep  # noqa: E402
from test_dsee_options import _default_state, _dsee_state, _track, _with_rgb  # noqa: E402


def _check(state):
    t = _track(state, mutate=_with_rgb)
    coord, n_prim, krow = eep_coordinate(t, state)
    prim = [int(x) for x in t.eep if x > 0]
    assert n_prim == len(prim) >= 4
    # exactly k_j at every primary, non-decreasing in between, NaN after the last primary
    assert np.allclose(coord[np.array(prim) - 1], krow[:len(prim)], atol=1e-12)
    span = coord[: prim[-1]]
    assert np.all(np.isfinite(span)) and np.all(np.diff(span) >= -1e-12)
    assert np.all(np.isnan(coord[prim[-1]:]))
    # reading the history at integer coordinates reproduces the .eep rows
    s = secondary_eep(copy.deepcopy(t), state)
    rows = np.arange(1, krow[len(prim) - 1] + 1)
    for name in ('log_L', 'log_Teff', 'star_age'):
        i = [c.name for c in t.cols].index(name)
        want = s.tr[i, :len(rows)]
        got = np.interp(rows, span, t.tr[i, : prim[-1]])
        assert np.allclose(got, want, rtol=1e-10, atol=1e-10), (name, np.max(np.abs(got - want)))
    return krow[:len(prim)]


def test_coord_default_rules():
    assert _check(_default_state())[:5] == [1, 202, 353, 454, 605]


def test_coord_dsee_rules():
    assert _check(_dsee_state())[:6] == [1, 202, 353, 454, 491, 605]


if __name__ == '__main__':
    test_coord_default_rules()
    test_coord_dsee_rules()
    print("OK")
