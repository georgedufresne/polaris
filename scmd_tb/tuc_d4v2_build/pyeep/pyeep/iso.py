"""Isochrone construction from EEP-based tracks.

Python port of ``src/make_iso.f90`` (plus the pieces of
``src/iso_eep_support.f90`` and MESA's ``interp_1d`` it depends on) from the
Fortran ``iso`` project.  This is Section 3.2 of Dotter (2016): for each age,
loop over the EEPs, find the tracks that define each EEP, interpolate in the
age--initial-mass relation to get the initial mass at that age, then
interpolate every other quantity in initial mass.

The Fortran loops over EEPs one at a time and, for every column, rebuilds a
full cubic interpolant across all tracks.  Here every step is vectorised
across *all EEPs of an isochrone at once*, and the interpolant is evaluated
locally (only the four tracks bracketing the target are touched).  The
floating-point operations are the same, in the same order, so the output
reproduces the Fortran bit-for-bit (see the module README for the validation).

The grid is read once and reused across ages: :class:`IsochroneGrid` holds the
tracks in padded NumPy arrays, and :meth:`IsochroneGrid.isochrone` builds one
age.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass, field
from typing import List, Optional, Sequence

import numpy as np

from .support import (
    STAR_HIGH_MASS,
    STAR_LABEL,
    STAR_LOW_MASS,
    UNKNOWN,
    AGE_SCALE_LINEAR,
    AGE_SCALE_LOG10,
    _fmt_1p_e13_5,
    _fmt_es32_16e3,
)

#: matches ``tiny`` in do_isochrone_for_age (mass matching tolerance)
MASS_TINY = 1.0e-12
#: matches ``tiny`` in num_binary_search.inc
SEARCH_TINY = 1.0e-13


class MakeIsoError(RuntimeError):
    """A condition on which the Fortran ``make_iso`` stops or writes nothing."""


# ---------------------------------------------------------------------------
# controls (the &iso_controls namelist)
# ---------------------------------------------------------------------------

@dataclass
class IsoControls:
    """``&iso_controls`` namelist parameters; defaults match make_iso.f90."""

    iso_debug: bool = False
    do_smooth: bool = True
    do_PAV: bool = True
    do_linear_interpolation: bool = False
    log_age_delta: float = 1.0
    very_low_mass_limit: float = 0.5
    set_max_eep_number: bool = False
    new_max_eep_number: int = -1

    @classmethod
    def from_namelist(cls, nml: dict) -> 'IsoControls':
        """Build from a lower-cased ``key -> value`` namelist dict."""
        c = cls()
        c.iso_debug = bool(nml.get('iso_debug', c.iso_debug))
        c.do_smooth = bool(nml.get('do_smooth', c.do_smooth))
        c.do_PAV = bool(nml.get('do_pav', c.do_PAV))
        c.do_linear_interpolation = bool(
            nml.get('do_linear_interpolation', c.do_linear_interpolation))
        c.log_age_delta = float(nml.get('log_age_delta', c.log_age_delta))
        c.very_low_mass_limit = float(
            nml.get('very_low_mass_limit', c.very_low_mass_limit))
        c.set_max_eep_number = bool(
            nml.get('set_max_eep_number', c.set_max_eep_number))
        c.new_max_eep_number = int(
            nml.get('new_max_eep_number', c.new_max_eep_number))
        return c


# ---------------------------------------------------------------------------
# reading an .eep file  ->  read_eep
# ---------------------------------------------------------------------------

@dataclass
class EEPTrack:
    """One ``.eep`` file, as read by make_iso."""

    filename: str = ''
    version_string: str = ''
    MESA_revision_number: int = 0
    initial_Y: float = 0.0
    initial_Z: float = 0.0
    Fe_div_H: float = 0.0
    alpha_div_Fe: float = 0.0
    v_div_vcrit: float = 0.0
    initial_mass: float = 0.0
    ntrack: int = 0
    neep: int = 0
    ncol: int = 0
    has_phase: bool = False
    star_type: int = UNKNOWN
    eep: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=int))
    col_names: List[str] = field(default_factory=list)
    #: shape (ncol, ntrack); rows before eep[0] are NaN (never read)
    tr: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    phase: Optional[np.ndarray] = None


def _fixed(line: str, a: int, b: int) -> str:
    return line[a:b] if len(line) >= b else line[a:].ljust(b - a)


def read_eep(path: str) -> EEPTrack:
    """Read an ``.eep`` track (header fields are fixed-width, as in Fortran)."""
    with open(path, 'r') as fh:
        lines = fh.read().splitlines()

    x = EEPTrack(filename=path)
    x.version_string = _fixed(lines[0], 25, 33)              # (25x,a8)
    x.MESA_revision_number = int(_fixed(lines[1], 25, 33))   # (25x,i8)
    ln = lines[4]                            # (2x,f6.4,1p1e13.5,0p3f9.2)
    x.initial_Y = float(_fixed(ln, 2, 8))
    x.initial_Z = float(_fixed(ln, 8, 21))
    x.Fe_div_H = float(_fixed(ln, 21, 30))
    x.alpha_div_Fe = float(_fixed(ln, 30, 39))
    x.v_div_vcrit = float(_fixed(ln, 39, 48))
    ln = lines[7]                            # (2x,1p1e16.10,3i8,a8,2x,a10)
    x.initial_mass = float(_fixed(ln, 2, 18))
    x.ntrack = int(_fixed(ln, 18, 26))
    x.neep = int(_fixed(ln, 26, 34))
    ncol = int(_fixed(ln, 34, 42))
    phase_info = _fixed(ln, 42, 50)
    type_label = _fixed(ln, 52, 62)
    for i, label in enumerate(STAR_LABEL):
        if type_label == label:
            x.star_type = i + 1

    x.has_phase = 'YES' in phase_info
    x.ncol = ncol - 1 if x.has_phase else ncol

    x.eep = np.array(lines[8][8:].split()[:x.neep], dtype=int)
    if len(x.eep) != x.neep:
        raise ValueError(f"{path}: expected {x.neep} EEPs, found {len(x.eep)}")

    x.col_names = lines[11][1:].split()
    first = int(x.eep[0])
    nrow = x.ntrack - first + 1
    width = ncol
    body = lines[12:12 + nrow]
    vals = np.array(' '.join(body).split(), dtype=np.float64)
    if vals.size != nrow * width:
        raise ValueError(f"{path}: expected {nrow} rows of {width} values, "
                         f"got {vals.size} values")
    vals = vals.reshape(nrow, width)

    x.tr = np.full((x.ncol, x.ntrack), np.nan)
    x.tr[:, first - 1:] = vals[:, :x.ncol].T
    if x.has_phase:
        x.phase = np.full(x.ntrack, np.nan)
        x.phase[first - 1:] = vals[:, x.ncol]
    return x


# ---------------------------------------------------------------------------
# small numerical kernels
# ---------------------------------------------------------------------------

def _log10(a: np.ndarray) -> np.ndarray:
    """Elementwise libm ``log10``, as gfortran calls it.

    NumPy's vectorised log10 differs from libm by 1 ULP on a few percent of
    inputs, which is enough to move an EEP across the age window.  The arrays
    involved are small (tracks x EEPs), so the scalar loop is cheap.
    """
    out = np.empty(a.shape)
    flat_in = a.ravel()
    flat_out = out.ravel()
    for i, v in enumerate(flat_in.tolist()):
        if v > 0.0:
            flat_out[i] = math.log10(v)
        elif v == 0.0:
            flat_out[i] = -math.inf
        else:
            flat_out[i] = math.nan
    return out


def pav(y: np.ndarray) -> np.ndarray:
    """Pool-adjacent-violators, a port of ``PAV`` in iso_eep_support.f90.

    Returns a new array; the pooled mean is summed left to right, as gfortran
    does for the ``sum`` intrinsic, so pooled values match bit-for-bit.
    """
    y = np.array(y, dtype=np.float64)
    n = y.size
    if n < 2:
        return y
    lo = np.arange(n)
    hi = np.arange(n)
    while True:
        d = y[1:] - y[:-1]
        neg = np.flatnonzero(d < 0.0)
        if neg.size == 0:
            break
        i = int(neg[0])
        start, last = int(lo[i]), int(hi[i + 1])
        acc = 0.0
        for v in y[start:last + 1].tolist():
            acc += v
        y[start:last + 1] = acc / float(last - start + 1)
        lo[start:last + 1] = start
        hi[start:last + 1] = last
    return y


def _binary_search(V: np.ndarray, cnt: np.ndarray, val: float) -> np.ndarray:
    """Vectorised ``binary_search`` (num_binary_search.inc) with hint = 1.

    ``V`` is (E, W): row e holds a 1-based vector of length ``cnt[e]``
    (padding beyond is ignored).  Returns the 1-based ``loc`` for every row,
    following exactly the Fortran's early returns and bisection path, so the
    answer agrees even where the vector is not monotonic.
    """
    E = V.shape[0]
    rows = np.arange(E)
    n = cnt
    v1 = V[:, 0]
    vn = V[rows, n - 1]
    v2 = V[rows, np.minimum(1, n - 1)]
    dec = vn < v1

    loc = np.full(E, -1, dtype=np.int64)
    done = n <= 1
    loc[done] = n[done]

    def settle(mask, value):
        nonlocal done
        m = mask & ~done
        loc[m] = value[m] if isinstance(value, np.ndarray) else value
        done = done | m

    settle(dec & (val > v1), 0)
    settle(~dec & (val < v1), 0)
    settle(np.abs(val - v1) < SEARCH_TINY, 1)
    settle(dec & (val <= vn), n)
    settle(~dec & (val >= vn), n)
    hint_ok = n > 1
    settle(hint_ok & dec & (v1 >= val) & (val > v2), 1)
    settle(hint_ok & ~dec & (v1 <= val) & (val < v2), 1)

    first = np.ones(E, dtype=np.int64)
    last = n - 1
    active = ~done
    while active.any():
        mid = (first + last) // 2
        midc = np.clip(mid, 1, V.shape[1] - 1)
        vm = V[rows, midc - 1]
        vm1 = V[rows, midc]
        go_up = np.where(dec, vm >= val, vm <= val)
        found = go_up & np.where(dec, val > vm1, val < vm1)
        hit = active & found
        loc[hit] = mid[hit]
        active = active & ~found
        step_up = active & go_up
        step_dn = active & ~go_up
        first = np.where(step_up, mid + 1, first)
        last = np.where(step_dn, mid - 1, last)
        active = active & (first <= last)
    return loc


def _interp_locate(X: np.ndarray, cnt: np.ndarray, xval: np.ndarray):
    """Interval search of MESA ``do_interp_values``.

    Returns the 1-based interval ``k`` and the offset ``delta = x - X[k]`` with
    ``x`` clipped to the range, scanning from the first interval exactly as the
    Fortran does.
    """
    E = X.shape[0]
    rows = np.arange(E)
    x1 = X[:, 0]
    xn = X[rows, cnt - 1]
    inc = x1 < xn
    xv = np.array(xval, dtype=np.float64, copy=True)
    xv = np.where(inc & (xv > xn), xn, np.where(inc & (xv < x1), x1, xv))
    xv = np.where(~inc & (xv < xn), xn, np.where(~inc & (xv > x1), x1, xv))

    k = np.ones(E, dtype=np.int64)
    xk = x1.copy()
    xkp1 = X[:, 1]
    active = np.ones(E, dtype=bool)
    while True:
        adv = active & ((inc & (xv > xkp1)) | (~inc & (xv < xkp1)))
        if not adv.any():
            break
        k = np.where(adv, k + 1, k)
        over = adv & (k >= cnt)
        if over.any():
            k = np.where(over, k - 1, k)
            xv = np.where(over, xkp1, xv)
            active = active & ~over
        move = adv & ~over
        xk = np.where(move, xkp1, xk)
        nxt = np.minimum(k, X.shape[1] - 1)  # 0-based index of point k+1
        xkp1 = np.where(move, X[rows, nxt], xkp1)
        active = active & move
    delta = xv - xk
    return k, delta


def _monotonic_mass_range(ages: np.ndarray, k: int):
    """0-based inclusive [lo, hi] of the decreasing run around crossing ``k``.

    Port of monotonic_mass_range (make_iso.f90), translated from 1-based.
    """
    n = ages.size
    lo = 0
    for j in range(k, -1, -1):
        if ages[j + 1] < ages[j]:
            lo = j
        else:
            break
    hi = n - 1
    for j in range(lo + 1, n):
        if ages[j - 1] > ages[j]:
            hi = j
        else:
            break
    return lo, hi


def _warn_repeated_x(eep: int, what: str) -> None:
    """Report an EEP dropped because the interpolation abscissae repeat.

    MESA's mk_pmcub refuses repeated x values (e.g. several truncated
    low-mass tracks that all end at exactly the same age).  The Fortran
    make_iso ignores that error and evaluates coefficients left in memory by
    the previous EEP; the port drops the point instead.
    """
    print(f" make_iso: EEP {eep} dropped: repeated x values in the "
          f"{what} interpolation")


def _sign(a):
    return np.copysign(1.0, a)


def _interior_slope(s_a, s_b, h_a, h_b):
    """Steffen slope at an interior point from neighbouring slopes/widths."""
    p = (s_a * h_b + s_b * h_a) / (h_a + h_b)
    m = np.minimum(np.minimum(np.abs(s_a), np.abs(s_b)), 0.5 * np.abs(p))
    return (_sign(s_a) + _sign(s_b)) * m


def _end_slope(s_near, s_far, h_near, h_far):
    """Steffen slope at an end point (mk_pmcub's p(1) / p(nx) branch)."""
    p = s_near * (1 + h_near / (h_near + h_far)) - s_far * h_near / (h_near + h_far)
    out = np.where(np.abs(p) > 2 * np.abs(s_near), 2 * s_near, p)
    return np.where(p * s_near <= 0, 0.0, out)


def steffen_eval(X: np.ndarray, cnt: np.ndarray, xval: np.ndarray,
                 fetch) -> np.ndarray:
    """Evaluate MESA's piecewise-monotonic (Steffen 1990) cubic, vectorised.

    Port of ``interp_pm`` + ``interp_value`` (interp_1d_pm.f90 /
    interp_1d_misc.f90).  Row ``e`` of ``X`` holds ``cnt[e]`` abscissae; the
    interpolant is evaluated at ``xval[e]``.  ``fetch(pos)`` must return the
    ordinates at 0-based positions ``pos`` (shape (E,)), as an array of shape
    (E,) or (E, m) -- this lets callers gather only the four tracks the local
    cubic needs instead of materialising every column of every track.

    Only the coefficients of the interval containing ``xval`` are formed,
    using the same expressions in the same order as mk_pmcub, so values agree
    with the full-array Fortran to the last bit.  Rows where two abscissae
    coincide get NaN (MESA flags this as an error; the Fortran caller then
    evaluates uninitialised coefficients).
    """
    E, W = X.shape
    rows = np.arange(E)
    k, delta = _interp_locate(X, cnt, xval)
    k0 = k - 1  # 0-based index of point k

    def pos(off):
        return np.clip(k0 + off, 0, np.maximum(cnt - 1, 0))

    xm1, x0, x1, x2 = (X[rows, pos(o)] for o in (-1, 0, 1, 2))
    fm1, f0, f1, f2 = (fetch(pos(o)) for o in (-1, 0, 1, 2))
    if f0.ndim == 2:
        col = (slice(None), None)
        xm1, x0, x1, x2, delta_b = (a[col] for a in (xm1, x0, x1, x2, delta))
        kk, nn = k[col], cnt[col]
    else:
        delta_b, kk, nn = delta, k, cnt

    with np.errstate(divide='ignore', invalid='ignore'):
        hm1 = x0 - xm1
        h0 = x1 - x0
        h1 = x2 - x1
        sm1 = (f0 - fm1) / hm1
        s0 = (f1 - f0) / h0
        s1 = (f2 - f1) / h1

        # slope at point k: start-end formula if k == 1, else interior
        d_k = np.where(kk == 1, _end_slope(s0, s1, h0, h1),
                       _interior_slope(sm1, s0, hm1, h0))
        # slope at point k+1: finish-end formula if k+1 == n, else interior
        d_k1 = np.where(kk + 1 == nn, _end_slope(s0, sm1, h0, hm1),
                        _interior_slope(s0, s1, h0, h1))

        c3 = (3 * s0 - 2 * d_k - d_k1) / h0
        c4 = (d_k + d_k1 - 2 * s0) / (h0 * h0)
        cubic = f0 + delta_b * (d_k + delta_b * (c3 + delta_b * c4))
        # nx == 2: mk_pmlinear sets f(2,1)=s, f(3,:)=f(4,:)=0
        zero = np.zeros_like(s0)
        linear = f0 + delta_b * (s0 + delta_b * (zero + delta_b * zero))
    out = np.where(nn == 2, linear, cubic)

    # mk_pmcub rejects any repeated abscissa anywhere in the vector
    j = np.arange(W - 1)
    within = j[None, :] < (cnt - 1)[:, None]
    dup = ((X[:, 1:] - X[:, :-1]) == 0) & within
    bad = dup.any(axis=1)
    if bad.any():
        out = np.array(out, dtype=np.float64, copy=True)
        out[bad] = np.nan
    return out


def _smooth_ages(M: np.ndarray, A: np.ndarray, cnt: np.ndarray,
                 rows_ok: np.ndarray) -> np.ndarray:
    """Vectorised ``smooth(masses, ages)``: inverse-distance weighted means.

    Applied in place in the Fortran's order -- point 2, point n-1, then points
    3..n-2 left to right -- so each step sees already-smoothed neighbours.
    Rows with fewer than 8 points are left alone, as in the Fortran.
    """
    A = A.copy()
    r_all = np.flatnonzero(rows_ok & (cnt >= 8))
    if r_all.size == 0:
        return A
    n_all = cnt[r_all]

    def npoint(r, idx_list, centre):
        x0 = M[r, centre]
        num = den = None
        for i in idx_list:                      # sum(y*w)/sum(w), left to right
            w = 1.0 / (1.0e-1 + np.abs(M[r, i] - x0))
            yw = A[r, i] * w
            num = yw if num is None else num + yw
            den = w if den is None else den + w
        return num / den

    one = np.ones_like(n_all)
    A[r_all, 1] = npoint(r_all, [one - 1, one, one + 1], one)
    A[r_all, n_all - 2] = npoint(r_all, [n_all - 3, n_all - 2, n_all - 1], n_all - 2)
    for i in range(2, int(n_all.max()) - 2):    # 0-based centres 2..n-3
        r = r_all[i <= n_all - 3]
        if r.size == 0:
            break
        A[r, i] = npoint(r, [i - 2, i - 1, i, i + 1, i + 2], i)
    return A


# ---------------------------------------------------------------------------
# the isochrone grid
# ---------------------------------------------------------------------------

@dataclass
class Isochrone:
    """One isochrone: EEP numbers, quantities per EEP and (optional) phase."""

    age: float
    age_scale: int
    eep: np.ndarray            # (neep,) int
    data: np.ndarray           # (neep, ncol); column i_Minit is initial_mass
    col_names: List[str]
    phase: Optional[np.ndarray] = None
    initial_Y: float = 0.0
    initial_Z: float = 0.0
    Fe_div_H: float = 0.0
    alpha_div_Fe: float = 0.0
    v_div_vcrit: float = 0.0

    @property
    def neep(self) -> int:
        return int(self.eep.size)

    def column(self, name: str) -> np.ndarray:
        return self.data[:, self.col_names.index(name)]


class IsochroneGrid:
    """A set of EEP tracks of one composition, ready to make isochrones.

    ``tracks`` must be the good (non-ignored) tracks in increasing mass order;
    ``col_names`` is the history column list (``star_age`` included, it is
    renamed ``initial_mass`` in the output, as in the Fortran).
    """

    def __init__(self, tracks: Sequence[EEPTrack], col_names: Sequence[str],
                 controls: Optional[IsoControls] = None):
        if len(tracks) < 1:
            raise MakeIsoError('no valid EEP tracks')
        self.controls = controls or IsoControls()
        self.tracks = list(tracks)
        names = list(col_names)
        if 'star_age' not in names:
            raise MakeIsoError("history column list has no 'star_age'")
        self.i_age = names.index('star_age')
        self.i_Minit = self.i_age
        names[self.i_Minit] = 'initial_mass'
        self.col_names = names
        self.ncol = len(names)
        for t in self.tracks:
            if t.ncol != self.ncol:
                raise MakeIsoError(
                    f"{t.filename}: {t.ncol} data columns but the history "
                    f"column list has {self.ncol}")
        self.has_phase = self.tracks[0].has_phase
        # composition written to the output header; the driver overrides
        # these with the values from the input file, as make_iso does
        t0 = self.tracks[0]
        self.initial_Y, self.initial_Z = t0.initial_Y, t0.initial_Z
        self.Fe_div_H, self.alpha_div_Fe = t0.Fe_div_H, t0.alpha_div_Fe
        self.v_div_vcrit = t0.v_div_vcrit
        self._build()

    # -- set-up shared by every age ----------------------------------------
    def _build(self) -> None:
        c = self.controls
        s = self.tracks
        n = len(s)
        self.n = n
        max_eep = max(t.ntrack for t in s)
        if (c.set_max_eep_number and c.new_max_eep_number > 0
                and c.new_max_eep_number <= max_eep):
            max_eep = c.new_max_eep_number
        self.max_eep = max_eep

        self.M = np.array([t.initial_mass for t in s], dtype=np.float64)
        self.eep_first = np.array([int(t.eep[0]) for t in s])
        self.eep_last = np.array([int(t.eep[t.neep - 1]) for t in s])

        # tracks incomplete for their type are dropped from every EEP
        max_neep_low = 0
        max_neep_high = 0
        for t in s:
            if t.star_type < STAR_HIGH_MASS:
                cur = max_neep_low
            else:
                cur = max_neep_high
            if c.set_max_eep_number and t.eep[t.neep - 1] > c.new_max_eep_number:
                hits = np.flatnonzero(t.eep == c.new_max_eep_number)
                if hits.size:
                    cur = int(hits[0]) + 1
            else:
                cur = max(cur, t.neep)
            if t.star_type < STAR_HIGH_MASS:
                max_neep_low = cur
            else:
                max_neep_high = cur
        base = np.zeros(n, dtype=bool)
        for k, t in enumerate(s):
            if t.star_type == STAR_HIGH_MASS and t.neep < max_neep_high:
                base[k] = True
            elif (t.initial_mass > c.very_low_mass_limit
                  and t.star_type == STAR_LOW_MASS and t.neep < max_neep_low):
                base[k] = True
        self.skip_base = base

        # star_age and its log at every (track, EEP); NaN where undefined
        E = max_eep
        age = np.full((n, E), np.nan)
        for k, t in enumerate(s):
            m = min(t.ntrack, E)
            age[k, :m] = t.tr[self.i_age, :m]
        self.age = age
        self.log_age = _log10(age)
        eeps = np.arange(1, E + 1)
        self.defined = ((eeps[None, :] >= self.eep_first[:, None])
                        & (eeps[None, :] <= self.eep_last[:, None]))

        # all track data, (n, E, ncol), for gathering
        TR = np.full((n, E, self.ncol), np.nan)
        for k, t in enumerate(s):
            m = min(t.ntrack, E)
            TR[k, :m, :] = t.tr[:, :m].T
        self.TR = TR
        if self.has_phase:
            self.phase_tbl = [t.phase for t in s]

    # -- one age -----------------------------------------------------------
    def isochrone(self, age_value: float,
                  age_scale: int = AGE_SCALE_LINEAR,
                  double_eep: bool = False) -> Isochrone:
        """Build the isochrone for one age (``do_isochrone_for_age``)."""
        age = math.log10(age_value) if age_scale == AGE_SCALE_LINEAR \
            else float(age_value)
        if double_eep:
            result, result2, valid = self._double(age)
        else:
            result, valid = self._single(age)
            result2 = None
        return self._assemble(age_value, age_scale, result, result2, valid)

    def _single(self, age: float):
        """The standard method: at most one point per EEP."""
        c = self.controls
        n, E = self.n, self.max_eep
        max_age = age + c.log_age_delta
        min_age = age - c.log_age_delta

        with np.errstate(invalid='ignore'):
            skip = (self.skip_base[:, None] | ~self.defined
                    | (self.log_age > max_age) | (self.log_age < min_age))
        # the Fortran's top-down monotonicity filter: its inner loop
        # `do l=k-1,1` only runs for k=2, so it only ever compares the two
        # lowest-mass tracks.  Reproduced as written.
        if n >= 2:
            drop = (~skip[1]) & (~skip[0]) & (self.age[1] > self.age[0])
            skip[0] |= drop

        keep = ~skip
        count = keep.sum(axis=0)
        ecols = np.flatnonzero(count >= 2)           # 0-based EEP index
        result = np.zeros((E, self.ncol))
        valid = np.zeros(E, dtype=bool)

        if ecols.size:
            # compacted per-EEP vectors: kept tracks first, in mass order
            order = np.argsort(skip[:, ecols], axis=0, kind='stable').T  # (Ev, n)
            cnt = count[ecols]
            rows = np.arange(ecols.size)
            ages_raw = self.log_age[order, ecols[:, None]]
            masses = self.M[order]
            W = n
            pad = np.arange(W)[None, :] >= cnt[:, None]
            ages_raw = np.where(pad, ages_raw[rows, (cnt - 1)][:, None], ages_raw)
            masses = np.where(pad, masses[rows, (cnt - 1)][:, None], masses)

            if c.do_smooth:
                ok = np.all(np.where(pad, True, masses > c.very_low_mass_limit),
                            axis=1)
                ages1 = _smooth_ages(masses, ages_raw, cnt, ok)
                ages1 = np.where(pad, ages1[rows, (cnt - 1)][:, None], ages1)
            else:
                ages1 = ages_raw

            # ---- pass 1: initial mass from age -------------------------
            loc = _binary_search(ages1, cnt, age)
            inrange = (loc >= 1) & (loc <= cnt - 1)

            # monotonic(masses): with n > 2 the masses must not reverse
            dm = masses[:, 1:] - masses[:, :-1]
            live = np.arange(W - 1)[None, :] < (cnt - 1)[:, None]
            asc = masses[:, 0] <= masses[rows, cnt - 1]
            viol = np.where(asc[:, None], dm < 0, dm > 0) & live
            nonmono = viol.any(axis=1) & (cnt > 2) & inrange
            if nonmono.any():
                raise MakeIsoError(
                    f' masses not monotonic in do_isochrone_for_age: {age!r}')

            sel = np.flatnonzero(inrange)
            mass = np.zeros(ecols.size)
            ok1 = np.zeros(ecols.size, dtype=bool)
            if sel.size:
                L = loc[sel] - 1                       # 0-based
                if c.do_linear_interpolation:
                    x1 = ages1[sel, L]
                    x2 = ages1[sel, L + 1]
                    y1 = masses[sel, L]
                    y2 = masses[sel, L + 1]
                    dx = x2 - x1
                    mass[sel] = ((x2 - age) / dx) * y1 + ((age - x1) / dx) * y2
                    ok1[sel] = True
                else:
                    Ms = masses[sel]
                    val = steffen_eval(ages1[sel], cnt[sel],
                                       np.full(sel.size, age),
                                       lambda p: Ms[np.arange(sel.size), p])
                    good = ~np.isnan(val)
                    for e in ecols[sel[~good]].tolist():
                        _warn_repeated_x(e + 1, 'mass from age')
                    mass[sel[good]] = val[good]
                    ok1[sel[good]] = True

            vidx = np.flatnonzero(ok1)
            if c.do_PAV and vidx.size:
                mass[vidx] = pav(mass[vidx])
            result[ecols[vidx], self.i_Minit] = mass[vidx]
            valid[ecols[vidx]] = True

            # ---- pass 2: everything else from initial mass -------------
            if c.do_smooth:
                loc2 = _binary_search(ages_raw, cnt, age)
            else:
                loc2 = loc
            in2 = ok1 & (loc2 >= 1) & (loc2 <= cnt - 1)
            sel = np.flatnonzero(in2)
            if sel.size:
                others = [i for i in range(self.ncol) if i != self.i_Minit]
                ee = ecols[sel]
                mm = mass[sel]
                if c.do_linear_interpolation:
                    L = loc2[sel] - 1
                    m1 = masses[sel, L]
                    m2 = masses[sel, L + 1]
                    klo = self._match_track(m1)
                    khi = self._match_track(m2)
                    Y1 = self.TR[klo, ee][:, others]
                    Y2 = self.TR[khi, ee][:, others]
                    dx = m2 - m1
                    alfa = ((m2 - mm) / dx)[:, None]
                    beta = ((mm - m1) / dx)[:, None]
                    vals = alfa * Y1 + beta * Y2
                    result[ee[:, None], np.array(others)[None, :]] = vals
                else:
                    trk = order[sel]                      # (S, n)
                    TRo = self.TR
                    ar = np.arange(sel.size)

                    def fetch(p):
                        return TRo[trk[ar, p], ee][:, others]

                    vals = steffen_eval(masses[sel], cnt[sel], mm, fetch)
                    bad = np.isnan(vals).any(axis=1) & self._dup_mass(masses[sel], cnt[sel])
                    for e in ee[bad].tolist():
                        _warn_repeated_x(e + 1, 'initial mass')
                    vals[bad] = 0.0
                    valid[ee[bad]] = False
                    result[ee[:, None], np.array(others)[None, :]] = vals

        return result, valid.astype(int)

    def _double(self, age: float):
        """``eep_style = double``: up to two points per EEP (Dotter 2016, 3.4).

        A scalar per-EEP loop; this mode is rarely used.  It mirrors the
        Fortran, including its use of the compacted index ``k`` to look up
        ``skip(k, eep)`` -- a track index -- in the pass loops.
        """
        c = self.controls
        n, E = self.n, self.max_eep
        skip = self.skip_base[:, None] | ~self.defined
        result1 = np.zeros((E, self.ncol))
        result2 = np.zeros((E, self.ncol))
        valid = np.zeros(E, dtype=int)
        others = [i for i in range(self.ncol) if i != self.i_Minit]
        one = np.ones(1, dtype=np.int64)

        def crossings(a):
            return [k for k in range(len(a) - 1)
                    if not skip[k, e] and a[k] > a[k + 1]
                    and a[k] >= age and a[k + 1] < age]

        def locate(a):
            return int(_binary_search(a[None, :], one * len(a), age)[0])

        for e in range(E):
            idx = np.flatnonzero(~skip[:, e])
            cnt = idx.size
            if cnt < 2:
                continue
            ages_raw = self.log_age[idx, e]
            masses = self.M[idx]
            ages = ages_raw
            if c.do_smooth and np.all(masses > c.very_low_mass_limit):
                ages = _smooth_ages(masses[None, :], ages_raw[None, :],
                                    one * cnt, np.ones(1, dtype=bool))[0]
            loc = locate(ages)
            if loc < 1 or loc > cnt - 1:
                continue
            if cnt > 2:
                asc = masses[0] <= masses[-1]
                d = np.diff(masses)
                if (d < 0).any() if asc else (d > 0).any():
                    raise MakeIsoError(
                        f' masses not monotonic in do_isochrone_for_age: {age!r}')
            for p, k in enumerate(crossings(ages)[:2]):
                lo, hi = _monotonic_mass_range(ages, k)
                xa = ages[lo:hi + 1]
                m = steffen_eval(xa[None, :], one * xa.size, np.full(1, age),
                                 lambda q: masses[lo:hi + 1][q])[0]
                if np.isnan(m):
                    _warn_repeated_x(e + 1, 'mass from age')
                    break
                (result1 if p == 0 else result2)[e, self.i_Minit] = m
                valid[e] = p + 1

        for e in range(E):
            if valid[e] < 1:
                continue
            idx = np.flatnonzero(~skip[:, e])
            cnt = idx.size
            ages = self.log_age[idx, e]
            masses = self.M[idx]
            loc = locate(ages)
            if loc < 1 or loc > cnt - 1:
                continue
            for p, _ in enumerate(crossings(ages)[:2]):
                res = result1 if p == 0 else result2
                v = steffen_eval(masses[None, :], one * cnt,
                                 np.full(1, res[e, self.i_Minit]),
                                 lambda q: self.TR[idx[q], e][:, others])[0]
                if np.isnan(v).any() and self._dup_mass(masses[None, :], one * cnt)[0]:
                    _warn_repeated_x(e + 1, 'initial mass')
                    valid[e] = 0
                    break
                res[e, others] = v
                valid[e] = p + 1
        return result1, result2, valid

    def _assemble(self, age_value, age_scale, result1, result2, valid):
        rows, eeps = [], []
        for e in np.flatnonzero(valid > 0).tolist():
            rows.append(result1[e]); eeps.append(e + 1)
            if valid[e] == 2:
                rows.append(result2[e]); eeps.append(e + 1)
        data = np.array(rows).reshape(len(rows), self.ncol)
        eep_arr = np.array(eeps, dtype=int)
        iso = Isochrone(age=age_value, age_scale=age_scale, eep=eep_arr,
                        data=data, col_names=self.col_names,
                        initial_Y=self.initial_Y, initial_Z=self.initial_Z,
                        Fe_div_H=self.Fe_div_H, alpha_div_Fe=self.alpha_div_Fe,
                        v_div_vcrit=self.v_div_vcrit)
        if self.has_phase:
            iso.phase = self._phases(data[:, self.i_Minit], eep_arr)
        return iso

    # -- helpers -----------------------------------------------------------
    def _match_track(self, m: np.ndarray) -> np.ndarray:
        """Last track whose initial mass is within 1e-12 (the klo/khi loop)."""
        hit = np.abs(m[:, None] - self.M[None, :]) < MASS_TINY
        last = self.n - 1 - np.argmax(hit[:, ::-1], axis=1)
        return last

    @staticmethod
    def _dup_mass(X: np.ndarray, cnt: np.ndarray) -> np.ndarray:
        W = X.shape[1]
        within = np.arange(W - 1)[None, :] < (cnt - 1)[:, None]
        return (((X[:, 1:] - X[:, :-1]) == 0) & within).any(axis=1)

    def _phases(self, minit: np.ndarray, eeps: np.ndarray) -> np.ndarray:
        """Phase of the first track heavier than each point (or the last)."""
        out = np.empty(minit.size)
        Mlim = self.M[:self.n - 1]
        for j, (m, e) in enumerate(zip(minit.tolist(), eeps.tolist())):
            above = np.flatnonzero(m < Mlim)
            k = int(above[0]) if above.size else self.n - 1
            t = self.tracks[k]
            out[j] = t.phase[min(t.ntrack, e) - 1]
        return out


# ---------------------------------------------------------------------------
# writing  ->  write_isochrones_to_file
# ---------------------------------------------------------------------------

_DASHES = '# ' + '-' * 86


def _fmt_es32_16e3_full(x: float) -> str:
    """``1pes32.16e3`` including gfortran's spelling of non-finite values."""
    if math.isnan(x):
        return 'NaN'.rjust(32)
    if math.isinf(x):
        return ('Infinity' if x > 0 else '-Infinity').rjust(32)
    return _fmt_es32_16e3(x)


def format_rows(eeps: np.ndarray, values: np.ndarray) -> str:
    """Rows of ``(i5,299(1pes32.16e3))``, fast path through C formatting.

    ``%31.16E`` gives a two-digit exponent one character short of the
    Fortran's three; inserting the 0 restores the width.  Values that would
    need a three-digit exponent (or are non-finite) take the slow,
    per-value path.
    """
    if values.size == 0:
        return ''
    a = np.abs(values)
    with np.errstate(invalid='ignore'):
        risky = (~np.isfinite(values)) | (a >= 9.9e99) | ((a < 1.1e-99) & (a > 0))
    if risky.any():
        out = []
        for e, row in zip(eeps.tolist(), values.tolist()):
            out.append(f"{e:5d}" + ''.join(_fmt_es32_16e3_full(v) for v in row))
        return '\n'.join(out) + '\n'
    fmt = '%5d' + '%31.16E' * values.shape[1]
    text = '\n'.join(fmt % ((e,) + tuple(row))
                     for e, row in zip(eeps.tolist(), values.tolist())) + '\n'
    return text.replace('E+', 'E+0').replace('E-', 'E-0')


def isochrone_block(iso: Isochrone) -> str:
    """Text of one isochrone (write_isochrone_to_file_{orig,phase})."""
    has_phase = iso.phase is not None
    my_ncol = len(iso.col_names) + (3 if has_phase else 2)
    lines = [f"# number of EEPs, cols = {iso.neep:5d}{my_ncol:5d}",
             '#' + f"{1:4d}" + ''.join(f"{i:32d}" for i in range(2, my_ncol + 1))]
    age_name = ('log10_isochrone_age_yr' if iso.age_scale == AGE_SCALE_LOG10
                else 'isochrone_age_yr')
    names = [age_name] + list(iso.col_names) + (['phase'] if has_phase else [])
    lines.append('# EEP' + ''.join(f"{nm:>32}" for nm in names))
    cols = [np.full((iso.neep, 1), float(iso.age)), iso.data]
    if has_phase:
        cols.append(iso.phase[:, None])
    values = np.hstack(cols)
    return '\n'.join(lines) + '\n' + format_rows(iso.eep, values)


def write_isochrones(path: str, isos: Sequence[Isochrone],
                     version_string: str, MESA_revision_number: int,
                     initial_Y: float, initial_Z: float, Fe_div_H: float,
                     alpha_div_Fe: float, v_div_vcrit: float) -> None:
    """Write an isochrone set in the Fortran make_iso text format."""
    head = [
        f"# MIST version number  = {version_string:<8.8}",
        f"# MESA revision number = {MESA_revision_number:8d}",
        _DASHES,
        '#  Yinit        Zinit   [Fe/H]   [a/Fe]  v/vcrit'.ljust(88),
        f"# {initial_Y:6.4f}{_fmt_1p_e13_5(initial_Z)}"
        f"{Fe_div_H:9.2f}{alpha_div_Fe:9.2f}{v_div_vcrit:9.2f}",
        _DASHES,
        f"# number of isochrones = {len(isos):5d}",
        _DASHES,
    ]
    parts = ['\n'.join(head) + '\n']
    for i, iso in enumerate(isos):
        parts.append(isochrone_block(iso))
        if i < len(isos) - 1:
            parts.append('\n\n')
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    with open(path, 'w') as fh:
        fh.write(''.join(parts))
