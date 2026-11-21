"""FSPS phase assignment for EEP tracks.

Port of ``src/phases.f90``.  Each point in the EEP track is tagged with an
integer-valued FSPS "phase" derived from which primary-EEP interval it falls
in, with a special case promoting hot, H-poor high-mass stars to Wolf-Rayet.
"""

from __future__ import annotations

import numpy as np

from .support import EEPState, Track, STAR_HIGH_MASS

# FSPS phase codes
PRE_MS = -1.0
MS = 0.0
SGB = 1.0
RGB = 2.0
RED_HB_CLUMP = 3.0
EARLY_AGB = 4.0
TP_AGB = 5.0
POST_AGB = 6.0
BLUE_STRAGGLER = 7.0
BLUE_HB = 8.0
WOLF_RAYET = 9.0
UNDEFINED = -9.0


def set_track_phase(t: Track, state: EEPState) -> None:
    """Assign an FSPS phase to every point of the (secondary) EEP track ``t``."""
    n = t.ntrack
    phase = np.full(n, UNDEFINED)

    logTe = t.tr[state.i_logTe]
    have_surfH = state.i_surfH >= 0
    surfH = t.tr[state.i_surfH] if have_surfH else None

    for i in range(1, n + 1):  # 1-based
        # find nearest primary EEP index i_p such that i <= cumulative boundary
        next_eep = 0
        i_p = state.primary - 1  # value of j if the loop completes
        for j in range(1, state.primary):  # 1 .. primary-1
            next_eep += 1 + state.eep_interval[j - 1]
            if i <= next_eep:
                i_p = j
                break

        if state.add_RGBBase and i_p >= 5:
            i_p -= 1   # TAMS->RGBBase and RGBBase->RGBTip are both RGB
        if i_p == 1:
            ph = PRE_MS
        elif 2 <= i_p <= 3:
            ph = MS
        elif i_p == 4:
            ph = RGB
        elif 5 <= i_p <= 6:
            ph = RED_HB_CLUMP
        elif i_p == 7:
            ph = EARLY_AGB
        elif i_p == 8:
            ph = TP_AGB
        elif 9 <= i_p <= 10:
            ph = POST_AGB
        else:
            ph = UNDEFINED

        if (have_surfH and logTe[i - 1] >= 4.0 and surfH[i - 1] <= 3e-1 and
                t.star_type == STAR_HIGH_MASS):
            ph = WOLF_RAYET

        phase[i - 1] = ph

    t.phase = phase
    t.has_phase = True
