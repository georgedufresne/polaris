"""Fractional EEP coordinate for every model of a track.

``make_eep`` resamples a history file onto a fixed set of EEP rows. Emulator
training data instead keeps *every* MESA model, so each model needs its own
position on the EEP scale. This inverts the secondary-EEP placement:
``secondary_eep`` puts the secondaries between primaries ``j`` and ``j+1`` at
equal steps of the distance metric ``t.dist``, so a model at distance ``d``
between them has the fractional ``.eep`` row

    coord = k_j + (d - d_j) / (d_{j+1} - d_j) * (N_j + 1)

where ``k_j`` is the ``.eep`` row of primary ``j`` and ``N_j`` the number of
secondaries in that interval. The coordinate equals ``k_j`` at each primary, is
continuous and non-decreasing, and reading any history column at integer
coordinates reproduces the rows of the ``.eep`` track exactly (both are linear
in the distance). It depends on the same controls as the ``.eep`` files, so it
changes with ``input.nml`` / the EEP-interval file (see "Optional primary-EEP
rules" in the README).
"""

from __future__ import annotations

from typing import List, Tuple

import numpy as np

from .eep import primary_eep
from .support import EEPState, Track, distance_along_track


def eep_coordinate(t: Track, state: EEPState) -> Tuple[np.ndarray, int, List[int]]:
    """Return ``(coord, n_primary, primary_rows)`` for every model of ``t``.

    ``coord[i]`` is the fractional ``.eep`` row (1-based) of model ``i``
    (0-based array); models past the last primary EEP are NaN. ``n_primary`` is
    the number of primary EEPs found and ``primary_rows`` their ``.eep`` rows.
    Runs ``distance_along_track`` and ``primary_eep`` on ``t``. Raises
    ``ValueError`` if the track is too short or has fewer than two primaries.
    """
    if t.ntrack <= 3:
        raise ValueError(f"track too short for a distance metric (ntrack={t.ntrack})")
    distance_along_track(t, state)
    primary_eep(t, state)
    prim = [int(x) for x in t.eep if x > 0]
    if len(prim) < 2:
        raise ValueError("fewer than 2 primary EEPs found")

    d = np.asarray(t.dist, dtype=np.float64)
    coord = np.full(t.ntrack, np.nan, dtype=np.float64)
    krow = [1]
    for j in range(len(prim) - 1):
        krow.append(krow[-1] + int(state.eep_interval[j]) + 1)
    for j in range(len(prim) - 1):
        lo, hi = prim[j], prim[j + 1]
        seg = np.arange(lo - 1, hi)
        span = d[hi - 1] - d[lo - 1]
        if span <= 0.0:                      # identical distances: spread evenly over the models
            frac = np.linspace(0.0, 1.0, len(seg)) if len(seg) > 1 else np.zeros(1)
        else:
            frac = (d[seg] - d[lo - 1]) / span
        coord[seg] = krow[j] + frac * (int(state.eep_interval[j]) + 1)
    return coord, len(prim), krow
