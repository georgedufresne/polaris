"""Primary/secondary EEP identification and interpolation.

Port of ``src/eep.f90`` (only the parts used by ``make_eep``).

Primary EEPs are characteristic points in a stellar track (ZAMS, TAMS, RGB
tip, ...).  Once located, a fixed number of evenly-spaced *secondary* EEPs is
interpolated between each adjacent pair along the cumulative-distance metric,
producing a track with a uniform point structure suitable for interpolation.

All ``*guess*`` arguments and returned locations are 1-based, matching the
Fortran source, so the arithmetic in ``primary_eep`` stays identical.  The
underlying numpy arrays are indexed 0-based via a small ``_v`` helper.
"""

from __future__ import annotations

import numpy as np

from .support import (
    EEPState,
    Track,
    binary_search,
    SUB_STELLAR,
    STAR_LOW_MASS,
    STAR_HIGH_MASS,
)


def _v(arr: np.ndarray, i1: int) -> float:
    """1-based scalar access into a 1-D numpy array."""
    return float(arr[i1 - 1])


# ---------------------------------------------------------------------------
# primary EEPs
# ---------------------------------------------------------------------------

def primary_eep(t: Track, state: EEPState) -> None:
    """Locate the primary EEPs of ``t`` (fills ``t.eep``, 1-based, 0 = missing).

    Sequence for a normal H star::

        PreMS -> ZAMS -> TAMS(0.35) -> TAMS(1e-12) -> RGBTip -> ZAHB -> TAHB
        then, low-mass: TPAGB -> PostAGB -> WDCS
             high-mass: CarbonBurn

    With ``add_RGBBase`` an RGBBase primary sits between TAMS(1e-12) and
    RGBTip. The IAMS threshold is ``IAMS_Xc`` (0.35 by default).

    ``Xc_HeL_ladder`` replaces everything above the ZAMS with the ladder::

        PreMS -> ZAMS -> 11 x Xc/Xc(ZAMS) -> TAMS(1e-12) -> 5 x log(L_He/L)

    See :func:`ladder_eep`.
    """
    t.eep = np.zeros(state.primary, dtype=int)
    inc = 1
    ieep = 1  # 1-based cursor into t.eep

    if t.he_star:
        t.eep[ieep - 1] = ZAHB(t, state, 1)
        if check(t, ieep):
            return
        ieep += 1
        t.eep[ieep - 1] = TAHB(t, state, 1e-4, t.eep[ieep - 2] + inc)
        if check(t, ieep):
            return
        ieep += 1
        ieep = _finish(t, state, ieep, inc)
    else:
        t.eep[ieep - 1] = preMS_fudge(1)
        if check(t, ieep):
            return
        ieep += 1
        t.eep[ieep - 1] = ZAMS(t, state, t.eep[ieep - 2] + inc)
        if check(t, ieep):
            return
        ieep += 1
        if state.Xc_HeL_ladder:
            ladder_eep(t, state, ieep, inc)
            return
        t.eep[ieep - 1] = TAMS(t, state, state.IAMS_Xc, t.eep[ieep - 2] + inc)
        if check(t, ieep):
            return
        ieep += 1
        t.eep[ieep - 1] = TAMS(t, state, 1e-12, t.eep[ieep - 2] + inc)
        if check(t, ieep):
            return
        ieep += 1
        if state.add_RGBBase:
            # the tip is found first, then the base between TAMS and the tip
            tams = t.eep[ieep - 2]
            tip = RGBTip(t, state, tams + inc)
            t.eep[ieep - 1] = RGBBase(t, state, tams, tip) if tip > 0 else 0
            if check(t, ieep):
                return
            ieep += 1
            t.eep[ieep - 1] = tip
        else:
            t.eep[ieep - 1] = RGBTip(t, state, t.eep[ieep - 2] + inc)
        if check(t, ieep):
            return
        ieep += 1
        t.eep[ieep - 1] = ZAHB(t, state, t.eep[ieep - 2] + inc)
        if check(t, ieep):
            return
        ieep += 1
        t.eep[ieep - 1] = TAHB(t, state, 1e-4, t.eep[ieep - 2] + inc)
        if check(t, ieep):
            return
        ieep += 1
        _finish(t, state, ieep, inc)


def _finish(t: Track, state: EEPState, ieep: int, inc: int) -> int:
    """Append the terminal EEPs that depend on the star type."""
    if t.star_type == STAR_LOW_MASS:
        t.eep[ieep - 1] = TPAGB(t, state, t.eep[ieep - 2] + inc)
        if check(t, ieep):
            return ieep
        ieep += 1
        t.eep[ieep - 1] = PostAGB(t, state, t.eep[ieep - 2] + inc)
        if check(t, ieep):
            return ieep
        ieep += 1
        t.eep[ieep - 1] = WDCS(t, state, t.eep[ieep - 2] + inc)
    elif t.star_type == STAR_HIGH_MASS:
        t.eep[ieep - 1] = CarbonBurn(t, state, t.eep[ieep - 2] + inc)
    return ieep


def ladder_eep(t: Track, state: EEPState, ieep: int, inc: int) -> int:
    """Place the Xc/HeL ladder above an already-placed PreMS and ZAMS.

    Burning rungs rather than morphology: eleven fractions of the ZAMS central
    hydrogen across the main sequence, the TAMS, then five decades of helium
    luminosity up the giant branch.  Both ladders are monotone in their own
    quantity, so searching each rung from the previous one gives the same rows
    as searching every rung from the ZAMS (or the TAMS) and costs less.

    Returns the cursor after the last primary that was placed.  A track that
    stops early simply ends the chain, leaving the remaining primaries at 0.
    """
    Xc_zams = _v(t.tr[state.i_Xc], t.eep[ieep - 2])
    if Xc_zams <= 0.0:
        return ieep

    for frac in state.Xc_ladder:
        t.eep[ieep - 1] = Xc_below(t, state, frac * Xc_zams, t.eep[ieep - 2] + inc)
        if check(t, ieep):
            return ieep
        ieep += 1

    t.eep[ieep - 1] = TAMS(t, state, 1e-12, t.eep[ieep - 2] + inc)
    if check(t, ieep):
        return ieep
    ieep += 1

    for level in state.HeL_ladder:
        t.eep[ieep - 1] = HeL_rung(t, state, level, t.eep[ieep - 2] + inc)
        if check(t, ieep):
            return ieep
        ieep += 1
    return ieep


def check(t: Track, i: int) -> bool:
    """A primary EEP is "bad" (stop the chain) if it is 0 or the last point."""
    return t.eep[i - 1] == 0 or t.eep[i - 1] == t.ntrack


# ---------------------------------------------------------------------------
# individual primary-EEP finders
# ---------------------------------------------------------------------------

def preMS_fudge(guess: int) -> int:
    return guess


def ZAMS(t: Track, state: EEPState, guess: int) -> int:
    tr = t.tr
    Xc = tr[state.i_Xc]
    Xmax = _v(Xc, 1)
    Xmin = Xmax - 1.0e-3
    my_guess = guess if 1 <= guess <= t.ntrack else 1

    i = my_guess
    Xci = _v(Xc, i)
    while Xci > Xmin and i < t.ntrack:
        i += 1
        Xci = _v(Xc, i)
    ZAMS1 = i

    # 1-based argmax of log_g over [1, ZAMS1]
    logg_slice = tr[state.i_logg, 0:ZAMS1]
    ZAMS3 = int(np.argmax(logg_slice)) + 1

    result = ZAMS3
    if t.star_type == SUB_STELLAR:
        result = int(np.argmax(tr[state.i_Tc, :])) + 1
    return result


def TAMS(t: Track, state: EEPState, Xmin: float, guess: int) -> int:
    tr = t.tr
    Xc = tr[state.i_Xc]
    if guess < 1:
        my_guess = 1
    elif guess >= t.ntrack:
        return 0
    else:
        my_guess = guess
    for i in range(my_guess, t.ntrack + 1):
        if _v(Xc, i) < Xmin:
            return i
    # a very-low-mass star whose main sequence outlives the run gets the anchor
    # at its last row; t.vlm_fallback marks it. VLM_fallback_min_age is the age
    # the run must reach (1.5e10 in the Fortran).
    result = 0
    if (t.initial_mass <= state.very_low_mass_limit and
            _v(tr[state.i_age], t.ntrack) >= state.VLM_fallback_min_age):
        result = t.ntrack
        t.vlm_fallback = True
    if t.merger and state.merger_fallback:
        result = t.ntrack
    return result


def Xc_below(t: Track, state: EEPState, Xmin: float, guess: int) -> int:
    """First row at or after ``guess`` whose central H has fallen to ``Xmin``.

    A rung of the main-sequence ladder.  Unlike :func:`TAMS` there is no
    very-low-mass or merger fallback: a star that never burns down this far
    just has no such EEP, and the chain stops there.
    """
    if guess >= t.ntrack:
        return 0
    my_guess = max(guess, 1)
    Xc = t.tr[state.i_Xc]
    for i in range(my_guess, t.ntrack + 1):
        if _v(Xc, i) <= Xmin:
            return i
    return 0


def HeL_rung(t: Track, state: EEPState, level: float, guess: int) -> int:
    """First row at or after ``guess`` with log(L_He/L_total) >= ``level``.

    A rung of the giant-branch ladder.  Helium luminosity is what the core is
    actually doing, so these rungs stay in the same physical place at every
    mass, where the HR-diagram tip does not: a track that is truncated before
    it reaches a rung has no EEP there rather than a spurious one at its last
    row.
    """
    if guess >= t.ntrack or state.i_logLHe < 0:
        return 0
    my_guess = max(guess, 1)
    ratio = t.tr[state.i_logLHe] - t.tr[state.i_logL]
    for i in range(my_guess, t.ntrack + 1):
        if _v(ratio, i) >= level:
            return i
    return 0


def RGBBase(t: Track, state: EEPState, tams: int, tip: int) -> int:
    """Subgiant-to-giant knee between TAMS and the RGB tip (1-based rows).

    Along TAMS -> tip, with x = fraction of the drop to the coolest Teff and
    y = fraction of the rise to the tip log L, the knee is the row that
    maximises x - y: where the track turns from mostly horizontal (subgiant)
    to mostly vertical (giant branch).  Returns a row strictly between TAMS and
    the tip, or 0 if there is no room.
    """
    if tams < 1 or tip <= tams + 2:
        return 0
    tr = t.tr
    T = tr[state.i_logTe, tams - 1:tip]
    L = tr[state.i_logL, tams - 1:tip]
    x = (T[0] - T) / max(T[0] - T.min(), 1e-9)
    y = (L - L[0]) / max(L[-1] - L[0], 1e-9)
    k = int(np.argmax(x - y))
    k = min(max(k, 1), len(T) - 2)
    return tams + k


def RGBTip(t: Track, state: EEPState, guess: int) -> int:
    """RGB tip: max log L or min log Teff while central Y is still near its
    starting value.  By default the earlier of the two (Fortran); with
    ``RGBTip_take_later`` the later, so metal-poor 2-4 Msun stars that ignite He
    fainter than they crossed the Hertzsprung gap get their tip at ignition.
    With ``RGBTip_require_He_ignition`` a run that never reaches
    log(L_He/L_H) > He_ignition_log_ratio after ``guess`` has no tip (0).
    """
    tr = t.tr
    if guess < 1 or guess >= t.ntrack:
        return 0
    if state.RGBTip_require_He_ignition:
        if state.i_logLHe < 0 or state.i_logLH < 0:
            raise ValueError("RGBTip_require_He_ignition needs log_LHe and log_LH in the history file")
        ratio = (tr[state.i_logLHe, guess - 1:t.ntrack]
                 - tr[state.i_logLH, guess - 1:t.ntrack])
        if not np.any(ratio > state.He_ignition_log_ratio):
            return 0
    my_guess = guess
    Yc = tr[state.i_Yc]
    logL = tr[state.i_logL]
    logTe = tr[state.i_logTe]
    Ymin = _v(Yc, my_guess) - 1e-2
    Lmax = -99.0
    Tmin = 6.0
    RGBTip1 = 0
    RGBTip2 = 0
    for i in range(my_guess, t.ntrack + 1):
        if _v(Yc, i) > Ymin and _v(logL, i) > Lmax:
            Lmax = _v(logL, i)
            RGBTip1 = i
        if _v(Yc, i) > Ymin and _v(logTe, i) < Tmin:
            Tmin = _v(logTe, i)
            RGBTip2 = i
    if state.RGBTip_take_later:
        return max(RGBTip1, RGBTip2)
    return min(RGBTip1, RGBTip2)


def ZAHB(t: Track, state: EEPState, guess: int) -> int:
    tr = t.tr
    Yc = tr[state.i_Yc]
    logLHe = tr[state.i_logLHe]
    Tc = tr[state.i_Tc]
    if guess < 1 or guess > t.ntrack:
        my_guess = 1
    elif guess == t.ntrack:
        return 0
    else:
        my_guess = guess
    Ymin = _v(Yc, my_guess) - 3e-2
    LHemax = -99.0
    my_guess_2 = my_guess
    for i in range(my_guess, t.ntrack + 1):
        if _v(Yc, i) > Ymin and _v(logLHe, i) > LHemax:
            LHemax = _v(logLHe, i)
            my_guess_2 = i
    Tmin = 99.0
    result = 0
    for i in range(my_guess_2, t.ntrack + 1):
        if _v(Yc, i) > Ymin and _v(Tc, i) < Tmin:
            Tmin = _v(Tc, i)
            result = i
    if result == guess:
        result = 0
    return result


def TAHB(t: Track, state: EEPState, Ymin: float, guess: int) -> int:
    tr = t.tr
    Yc = tr[state.i_Yc]
    if guess < 1 or guess > t.ntrack:
        my_guess = 1
    elif guess == t.ntrack:
        return 0
    else:
        my_guess = guess
    for i in range(my_guess, t.ntrack + 1):
        if _v(Yc, i) < Ymin:
            return i
    return 0


def TPAGB(t: Track, state: EEPState, guess: int) -> int:
    tr = t.tr
    Ymin = 1e-6
    HeShellMin = 1e-1
    Yc = tr[state.i_Yc]
    he_core = tr[state.i_he_core]
    co_core = tr[state.i_co_core]
    if guess < 1 or guess > t.ntrack:
        my_guess = 1
    elif guess == t.ntrack:
        return 0
    else:
        my_guess = guess
    for i in range(my_guess, t.ntrack + 1):
        yc = _v(Yc, i)
        he_shell = _v(he_core, i) - _v(co_core, i)
        if yc < Ymin and he_shell < HeShellMin:
            return i
    return 0


def PostAGB(t: Track, state: EEPState, guess: int) -> int:
    tr = t.tr
    core_mass_frac_limit = 8e-1
    Tc = tr[state.i_Tc]
    co_core = tr[state.i_co_core]
    mass = tr[state.i_mass]
    if guess < 1 or guess > t.ntrack:
        my_guess = 1
    elif guess == t.ntrack:
        return 0
    else:
        my_guess = guess

    Tc_now = _v(Tc, my_guess)
    Tc_end = _v(Tc, t.ntrack)

    if Tc_now > Tc_end:  # has TP-AGB
        for i in range(my_guess, t.ntrack + 1):
            core_mass_frac = _v(co_core, i) / _v(mass, i)
            if core_mass_frac > core_mass_frac_limit:
                return i
        return 0
    else:  # has no TP-AGB
        return 0


def WDCS(t: Track, state: EEPState, guess: int) -> int:
    tr = t.tr
    gamma = tr[state.i_gamma]
    if guess < 1:
        my_guess = 1
    elif guess >= t.ntrack:
        return 0
    else:
        my_guess = guess
    for i in range(my_guess, t.ntrack + 1):
        if _v(gamma, i) >= state.center_gamma_limit:
            return i
    return 0


def CarbonBurn(t: Track, state: EEPState, guess: int) -> int:
    tr = t.tr
    limit_XY = 1e-8
    Xc = tr[state.i_Xc]
    Yc = tr[state.i_Yc]
    Cc = tr[state.i_Cc]
    if guess < 1:
        my_guess = 1
    elif guess >= t.ntrack:
        return 0
    else:
        my_guess = guess
    result = 0
    for i in range(my_guess, t.ntrack + 1):
        if (_v(Xc, i) < limit_XY and _v(Yc, i) < limit_XY and
                _v(Cc, i) < state.center_carbon_limit):
            result = i
            break
    return result


# ---------------------------------------------------------------------------
# secondary EEPs
# ---------------------------------------------------------------------------

def secondary_eep(t: Track, state: EEPState) -> Track:
    """Build an EEP track ``s`` from ``t`` with evenly-distributed secondaries.

    Requires ``primary_eep`` to have been run on ``t``.
    """
    s = Track()
    s.version_string = t.version_string
    s.initial_Y = t.initial_Y
    s.initial_Z = t.initial_Z
    s.Fe_div_H = t.Fe_div_H
    s.alpha_div_Fe = t.alpha_div_Fe
    s.v_div_vcrit = t.v_div_vcrit
    s.star_type = t.star_type

    # count present primaries and the secondaries that precede them
    num_p = 0
    num_s = 0
    for j in range(1, state.primary + 1):
        if t.eep[j - 1] > 0:
            num_p += 1
            if j > 1:
                num_s += state.eep_interval[j - 2]

    s.filename = t.filename
    s.initial_mass = t.initial_mass
    s.MESA_revision_number = t.MESA_revision_number
    s.ntrack = num_p + num_s
    s.neep = num_p
    s.ncol = t.ncol
    s.cols = list(t.cols)

    s.tr = np.zeros((s.ncol, s.ntrack))
    s.dist = np.zeros(s.ntrack)
    s.eep = np.zeros(s.neep, dtype=int)

    k = 1  # 1-based write cursor into s
    for j in range(1, num_p + 1):
        src = t.eep[j - 1]  # 1-based index into t
        s.tr[:, k - 1] = t.tr[:, src - 1]
        s.dist[k - 1] = t.dist[src - 1]
        s.eep[j - 1] = k
        if j < num_p:
            _eep_interpolate(t, state, j, k, s)
            k = k + state.eep_interval[j - 1] + 1

    return s


def _eep_interpolate(t: Track, state: EEPState, j: int, k: int, s: Track) -> None:
    """Fill the secondary EEPs between primaries ``j`` and ``j+1`` (linear).

    ``j`` is the primary index (1-based), ``k`` the 1-based location in ``s`` of
    primary ``j``.  Only the linear branch of the Fortran routine is ported (the
    cubic branch there is dead code, guarded by ``linear = .true.``).
    """
    if t.eep[j - 1] == t.eep[j]:  # eep(j) == eep(j+1)
        return

    num_eep = state.eep_interval[j - 1]
    d = t.dist  # 1-based via _v / index-1

    dist_lo = _v(d, t.eep[j - 1])
    dist_hi = _v(d, t.eep[j])
    total = dist_hi - dist_lo
    delta = total / float(num_eep + 1)

    j0 = int(t.eep[j - 1])
    n = t.ntrack
    for i in range(1, num_eep + 1):
        new_dist = dist_lo + delta * float(i)
        j0 = binary_search(n, d, j0, new_dist)
        j1 = j0 + 1
        alpha = (_v(d, j1) - new_dist) / (_v(d, j1) - _v(d, j0))
        beta = 1.0 - alpha
        s.tr[:, k + i - 1] = alpha * t.tr[:, j0 - 1] + beta * t.tr[:, j1 - 1]
        s.dist[k + i - 1] = new_dist
