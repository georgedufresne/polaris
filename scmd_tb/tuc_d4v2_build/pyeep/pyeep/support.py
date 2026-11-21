"""Core data structures and I/O for the EEP pipeline.

This is a Python port of the relevant subset of ``src/iso_eep_support.f90``
from the Fortran ``iso`` project.  It provides:

* the :class:`Column`, :class:`Track` and :class:`EEPState` data holders,
* parsing of a MESA ``history_columns.list`` file,
* reading of a MESA history file into a :class:`Track`,
* the cumulative "distance along track" metric used to place EEPs,
* stellar-type classification, and
* writing an ``.eep`` track file in the same textual format as the Fortran code.

Only the pieces required by ``make_eep`` are ported.  The binary (``.bin``)
caching and the isochrone-related routines are intentionally omitted.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# module-level constants (mirroring the Fortran parameters)
# ---------------------------------------------------------------------------

COL_WIDTH = 32
FILE_PATH = 256

# age scales for isochrones (kept for parity, unused by make_eep)
AGE_SCALE_LINEAR = 0
AGE_SCALE_LOG10 = 1

# stellar types for handling primary eeps
UNKNOWN = 1        # for initialization only
SUB_STELLAR = 2    # no fusion = brown dwarf
STAR_LOW_MASS = 3  # ends as a WD
STAR_HIGH_MASS = 4 # does not end as a WD

# right-justified in a field of width 10, exactly as the Fortran star_label
STAR_LABEL = ['   unknown', 'substellar', '  low-mass', ' high-mass']

CHECK_INITIAL_MASS = True
MASS_EPS = 1.0e-6
OLD_CORE_MASS_NAMES = False

# column type tags
COLUMN_INT = 0
COLUMN_DBL = 1

MAX_COL = 180


@dataclass
class Column:
    """A single history column: name, type tag and 1-based location."""

    name: str
    type: int = COLUMN_DBL
    loc: int = 0


@dataclass
class EEPState:
    """Holds the global configuration that the Fortran code kept in module scope.

    This bundles the history-column table, the EEP-interval specification and
    all of the tunable ``eep_controls`` namelist parameters so they can be
    passed explicitly instead of relying on global state.
    """

    # column table (built from history_columns.list)
    cols: List[Column] = field(default_factory=list)
    ncol: int = 0

    # zero-based indices of the "special" columns used for EEP finding.
    # (-1 means "not present")
    i_age: int = -1
    i_mass: int = -1
    i_logLH: int = -1
    i_logLHe: int = -1
    i_logTe: int = -1
    i_logL: int = -1
    i_logg: int = -1
    i_Tc: int = -1
    i_Rhoc: int = -1
    i_Xc: int = -1
    i_Yc: int = -1
    i_he_core: int = -1
    i_co_core: int = -1
    i_Cc: int = -1
    i_gamma: int = -1
    i_surfH: int = -1

    # EEP interval specification (from input.eep)
    primary: int = 0                       # number of primary EEPs
    eep_interval: List[int] = field(default_factory=list)  # len primary-1

    # eep_controls namelist parameters (defaults match iso_eep_support.f90)
    center_gamma_limit: float = 1.0e2
    center_carbon_limit: float = 1.0e-4
    log_center_T_limit: float = 9.0
    high_mass_limit: float = 1.0e1
    very_low_mass_limit: float = 0.5

    weight_center_rho_T_by_Xc: bool = False
    Teff_scale: float = 2.0
    logL_scale: float = 0.125
    age_scale: float = 0.05
    Rhoc_scale: float = 1.0
    Tc_scale: float = 1.0

    # optional primary-EEP rules (eep_controls). The defaults reproduce the
    # Fortran make_eep exactly; run/input_dsee.nml switches them all on.
    IAMS_Xc: float = 0.35                  # central X that defines the IAMS
    RGBTip_require_He_ignition: bool = False   # no tip unless He ignites after TAMS
    He_ignition_log_ratio: float = -3.0    # ignition = log(L_He/L_H) above this
    RGBTip_take_later: bool = False        # tip = later (not earlier) of max L / min Teff
    add_RGBBase: bool = False              # extra primary between TAMS and RGBTip
    VLM_fallback_min_age: float = 1.5e10   # very-low-mass TAMS fallback needs this age
    merger_fallback: bool = True           # mergers get TAMS at their last row

    # The Xc/HeL ladder (eep_controls: xc_hel_ladder).  A different primary set
    # altogether: PreMS, ZAMS, eleven rungs of central hydrogen, TAMS, then five
    # rungs of helium luminosity.  IAMS, RGBBase and RGBTip are not placed, so
    # the switches above have no effect while this is on.  19 primaries, and the
    # intervals come from run/input_dsee_ladder.eep.
    Xc_HeL_ladder: bool = False
    # rungs of Xc / Xc(ZAMS); strictly decreasing
    Xc_ladder: Tuple[float, ...] = (0.95, 0.90, 0.80, 0.70, 0.60, 0.50,
                                    0.40, 0.30, 0.20, 0.10, 0.01)
    # rungs of log(L_He / L_total) after the TAMS; strictly increasing
    HeL_ladder: Tuple[float, ...] = (-10.0, -8.0, -6.0, -4.0, -2.0)

    # file-format spec (from input.format); kept for reference
    head: int = 40
    main: int = 40
    xtra: int = 1


@dataclass
class Track:
    """An evolutionary track (history file or the derived EEP track)."""

    filename: str = ''
    version_string: str = 'unknown'
    cols: List[Column] = field(default_factory=list)
    has_phase: bool = False
    merger: bool = False
    vlm_fallback: bool = False             # TAMS was placed by the very-low-mass fallback
    write_distance: bool = False
    he_star: bool = False

    ncol: int = 0
    ntrack: int = 0
    neep: int = 0
    MESA_revision_number: int = 0
    star_type: int = UNKNOWN

    eep: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=int))

    initial_mass: float = 0.0
    initial_Y: float = 0.0
    Fe_div_H: float = 0.0
    initial_Z: float = 0.0
    v_div_vcrit: float = 0.0
    alpha_div_Fe: float = 0.0

    # tr has shape (ncol, ntrack) to mirror the Fortran t% tr(i, j)
    tr: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))
    dist: np.ndarray = field(default_factory=lambda: np.zeros(0))
    phase: np.ndarray = field(default_factory=lambda: np.zeros(0))


# ---------------------------------------------------------------------------
# binary search (port of num_binary_search.inc)
# ---------------------------------------------------------------------------

def binary_search(n: int, vec: np.ndarray, hint: int, val: float) -> int:
    """Port of the MESA-derived ``binary_search``.

    ``vec`` is a 1-based conceptual array; ``hint``/return value are 1-based to
    match the Fortran caller.  Only the increasing-values branch is exercised by
    ``make_eep`` (distances are monotonically increasing) but both branches are
    implemented for fidelity.
    """
    tiny = 1.0e-13

    def v(k: int) -> float:  # 1-based access into the numpy array
        return float(vec[k - 1])

    if n <= 1:
        return n

    if v(n) < v(1):  # decreasing values
        if val > v(1):
            return 0
        if abs(val - v(1)) < tiny:
            return 1
        if val <= v(n):
            return n
        if 0 < hint < n:
            if v(hint) >= val > v(hint + 1):
                return hint
        first, last = 1, n - 1
        while first <= last:
            mid = (first + last) // 2
            if v(mid) >= val:
                if val > v(mid + 1):
                    return mid
                first = mid + 1
            else:
                last = mid - 1
        return -1
    else:  # increasing values
        if val < v(1):
            return 0
        if abs(val - v(1)) < tiny:
            return 1
        if val >= v(n):
            return n
        if 0 < hint < n:
            if v(hint) <= val < v(hint + 1):
                return hint
        first, last = 1, n - 1
        while first <= last:
            mid = (first + last) // 2
            if v(mid) <= val:
                if val < v(mid + 1):
                    return mid
                first = mid + 1
            else:
                last = mid - 1
        return -1


# ---------------------------------------------------------------------------
# history_columns.list parsing (process_history_columns)
# ---------------------------------------------------------------------------

def process_history_columns(history_columns_list: str) -> List[Column]:
    """Parse a MESA ``history_columns.list`` file into a list of :class:`Column`.

    Blank lines and pure-comment lines are ignored.  A trailing ``! comment`` is
    stripped from a data line.  A column is flagged as integer-typed if its name
    contains ``number`` or begins with ``num_`` (matching the Fortran heuristic).
    Spaces inside a name become underscores.
    """
    cols: List[Column] = []
    with open(history_columns_list, 'r') as fh:
        for raw in fh:
            line = raw.replace('\t', ' ')
            if len(line.strip()) == 0:
                continue
            line = line.lstrip()
            i = line.find('!')
            # text of the (possible) column name, before any comment
            body = line if i < 0 else line[:i]
            body = body.strip()
            if len(body) == 0:
                # pure comment line
                continue

            is_int = False
            if 'number' in body:
                is_int = True
            elif body.startswith('num_'):
                is_int = True

            name = body.replace(' ', '_')
            cols.append(Column(name=name,
                               type=COLUMN_INT if is_int else COLUMN_DBL,
                               loc=len(cols) + 1))
    return cols


def locate_column(cols: List[Column], col_name: str) -> int:
    """Return the zero-based index of ``col_name`` in ``cols`` or -1 if absent."""
    for idx, c in enumerate(cols):
        if c.name == col_name:
            return idx
    return -1


def setup_columns(history_columns_list: str, state: EEPState) -> None:
    """Populate ``state.cols`` and the special-column indices (setup_columns)."""
    cols = process_history_columns(history_columns_list)
    state.cols = cols
    state.ncol = len(cols)

    state.i_age = locate_column(cols, 'star_age')
    state.i_mass = locate_column(cols, 'star_mass')
    if state.i_mass < 0:
        # could be a binary run
        state.i_mass = locate_column(cols, 'star_1_mass')
    state.i_logLH = locate_column(cols, 'log_LH')
    state.i_logLHe = locate_column(cols, 'log_LHe')
    state.i_logTe = locate_column(cols, 'log_Teff')
    state.i_logL = locate_column(cols, 'log_L')
    state.i_logg = locate_column(cols, 'log_g')
    state.i_Tc = locate_column(cols, 'log_center_T')
    state.i_Rhoc = locate_column(cols, 'log_center_Rho')
    state.i_Xc = locate_column(cols, 'center_h1')
    state.i_Yc = locate_column(cols, 'center_he4')
    state.i_Cc = locate_column(cols, 'center_c12')
    state.i_gamma = locate_column(cols, 'center_gamma')
    state.i_surfH = locate_column(cols, 'surface_h1')
    if OLD_CORE_MASS_NAMES:
        state.i_he_core = locate_column(cols, 'h1_boundary_mass')
        state.i_co_core = locate_column(cols, 'he4_boundary_mass')
    else:
        state.i_he_core = locate_column(cols, 'he_core_mass')
        state.i_co_core = locate_column(cols, 'c_core_mass')


# ---------------------------------------------------------------------------
# EEP interval spec (input.eep)  ->  set_eep_interval
# ---------------------------------------------------------------------------

def set_eep_interval(eep_input_file: str, state: EEPState) -> None:
    """Read ``input.eep`` and set ``state.primary`` and ``state.eep_interval``.

    File layout (see input.eep)::

        # comment
        <primary>                 number of primary EEPs
        <j> <interval>            repeated primary-1 times
    """
    with open(eep_input_file, 'r') as fh:
        lines = fh.readlines()
    fh_iter = iter(lines)
    next(fh_iter)  # comment line
    primary = int(next(fh_iter).split()[0])
    state.primary = primary
    state.eep_interval = [0] * (primary - 1)
    for _ in range(primary - 1):
        parts = next(fh_iter).split()
        j = int(parts[0])
        interval = int(parts[1])
        state.eep_interval[j - 1] = interval


def set_default_eep_interval(state: EEPState,
                             primary: int = 10,
                             n_secondary: int = 50) -> None:
    """Fallback used when input.eep is missing (matches the Fortran default)."""
    state.primary = primary
    state.eep_interval = [n_secondary] * (primary - 1)


# ---------------------------------------------------------------------------
# reading a MESA history file  ->  read_history_file
# ---------------------------------------------------------------------------

def _parse_revision(raw: str) -> int:
    """Parse a MESA ``version_number`` header value into an integer.

    Old MESA revisions are plain integers (e.g. ``12000``).  Newer releases use
    a quoted string such as ``"r22.05.1"``; in that case we extract the leading
    run of digits (``22``) so an integer field can still be written.  If nothing
    numeric is found, 0 is returned.
    """
    s = raw.strip().strip('"').strip("'")
    try:
        return int(float(s))
    except ValueError:
        pass
    digits = ''
    for ch in s:
        if ch.isdigit():
            digits += ch
        elif digits:
            break
    return int(digits) if digits else 0


def read_history_file(t: Track, state: EEPState, history_dir: str) -> None:
    """Read a MESA history file into ``t`` and derive distance / star type.

    The Fortran reader parses fixed-width columns; here we tokenize on
    whitespace, which is a robust equivalent for standard MESA history output
    and yields identical numeric values.  The mapping from the desired history
    columns (``state.cols``) to positions in the file is done by matching the
    column-name header row.
    """
    path = f"{history_dir}/{t.filename}"
    with open(path, 'r') as fh:
        lines = [ln.rstrip('\n') for ln in fh]

    # MESA history layout:
    #   0: header column numbers
    #   1: header names
    #   2: header values
    #   3: (blank)
    #   4: data column numbers
    #   5: data column names
    #   6...: data rows
    header_names = lines[1].split()
    header_values = lines[2].split()

    hdr = dict(zip(header_names, header_values))
    if 'version_number' in hdr:
        t.MESA_revision_number = _parse_revision(hdr['version_number'])
    if 'initial_mass' in hdr:
        t.initial_mass = float(hdr['initial_mass'])

    data_names = lines[5].split()
    data_rows = [ln for ln in lines[6:] if len(ln.strip()) > 0]

    ntrack = len(data_rows)
    ncol = state.ncol
    t.ncol = ncol
    t.ntrack = ntrack
    t.cols = list(state.cols)

    # position of each desired column within the file (or -1 -> filled with 0)
    file_index = {name: k for k, name in enumerate(data_names)}
    output = [file_index.get(c.name, -1) for c in state.cols]
    for c, o in zip(state.cols, output):
        if o < 0:
            print(f"do not have {c.name}")

    # parse data into tr[col, point]
    data = np.zeros((ntrack, len(data_names)))
    for j, row in enumerate(data_rows):
        vals = row.split()
        data[j, :len(vals)] = [float(x) for x in vals]

    tr = np.zeros((ncol, ntrack))
    for i, o in enumerate(output):
        if o >= 0:
            tr[i, :] = data[:, o]
    t.tr = tr

    t.dist = np.zeros(ntrack)
    distance_along_track(t, state)

    # The Fortran applies this correction *after* set_star_type_from_history,
    # but t% initial_mass has no default initializer there: when the history
    # header carries no initial_mass field it holds recycled heap memory.  We
    # resolve the mass from the track data first so the sub-stellar test is
    # deterministic and uses the real initial mass.
    if CHECK_INITIAL_MASS and state.i_mass >= 0:
        if abs(t.initial_mass - t.tr[state.i_mass, 0]) > MASS_EPS:
            t.initial_mass = float(t.tr[state.i_mass, 0])

    set_star_type_from_history(t, state)


# ---------------------------------------------------------------------------
# distance along track  ->  distance_along_track
# ---------------------------------------------------------------------------

def _sqdiff(x0: np.ndarray, x1: np.ndarray) -> np.ndarray:
    dx = x0 - x1
    return dx * dx


def distance_along_track(t: Track, state: EEPState) -> None:
    """Cumulative weighted distance used to distribute secondary EEPs."""
    tr = t.tr
    n = t.ntrack
    t.dist = np.zeros(n)
    if n <= 3:
        return

    if state.weight_center_rho_T_by_Xc:
        max_center_h1 = float(np.max(tr[state.i_Xc, :]))
        if max_center_h1 <= 0.0:
            max_center_h1 = 1.0
    else:
        max_center_h1 = 1.0

    logTe = tr[state.i_logTe]
    logL = tr[state.i_logL]
    Rhoc = tr[state.i_Rhoc]
    Tc = tr[state.i_Tc]
    # libm log10 (as used by gfortran), not NumPy's vectorised log10: the two
    # differ by 1 ULP on ~5% of inputs, which accumulates over the ~10^5-step
    # cumulative sum below and perturbs the interpolated secondary EEPs.
    # math.log10 raises on non-positive input where IEEE (and hence gfortran)
    # yields -inf / nan, so guard those out-of-domain ages explicitly.
    log_age = [math.log10(x) if x > 0.0 else (-math.inf if x == 0.0 else math.nan)
               for x in tr[state.i_age]]

    for j in range(1, n):
        if state.weight_center_rho_T_by_Xc:
            weight = max(0.0, tr[state.i_Xc, j] / max_center_h1)
        else:
            weight = 1.0

        tmp = state.Teff_scale * (logTe[j] - logTe[j - 1]) ** 2
        tmp += state.logL_scale * (logL[j] - logL[j - 1]) ** 2
        tmp += weight * state.Rhoc_scale * (Rhoc[j] - Rhoc[j - 1]) ** 2
        tmp += weight * state.Tc_scale * (Tc[j] - Tc[j - 1]) ** 2
        tmp += state.age_scale * (log_age[j] - log_age[j - 1]) ** 2

        t.dist[j] = t.dist[j - 1] + math.sqrt(tmp)


# ---------------------------------------------------------------------------
# stellar type classification  ->  set_star_type_from_history
# ---------------------------------------------------------------------------

def set_star_type_from_history(t: Track, state: EEPState) -> None:
    """Classify a track as sub-stellar / low-mass / high-mass (in that order)."""
    tr = t.tr
    n = t.ntrack  # 1-based count; last index is n-1

    # simple test for substellar: central H is essentially unchanged
    if (t.initial_mass < 0.1 and
            np.max(np.abs(tr[state.i_Xc, :] - tr[state.i_Xc, 0])) < 0.1):
        t.star_type = SUB_STELLAR
        return

    if state.i_gamma >= 0:
        # only reach center_gamma_limit if the star evolves to a WD
        if tr[state.i_gamma, n - 1] > state.center_gamma_limit:
            t.star_type = STAR_LOW_MASS
            return

    # simple test for high-mass stars: central C is depleted
    if (np.max(tr[state.i_Cc, :]) > 0.4 and
            tr[state.i_Cc, n - 1] < state.center_carbon_limit):
        t.star_type = STAR_HIGH_MASS
        return

    if state.i_Tc >= 0:
        if tr[state.i_Tc, n - 1] > state.log_center_T_limit:
            t.star_type = STAR_HIGH_MASS
        else:
            t.star_type = STAR_LOW_MASS

    # last-gasp test: the initial mass
    if t.initial_mass >= state.high_mass_limit:
        t.star_type = STAR_HIGH_MASS
    else:
        t.star_type = STAR_LOW_MASS


# ---------------------------------------------------------------------------
# writing an .eep track file  ->  write_track
# ---------------------------------------------------------------------------

def _fmt_es32_16e3(x: float) -> str:
    """Format a float like Fortran ``1pes32.16e3``.

    One nonzero digit before the decimal point, 16 after, a signed 3-digit
    exponent, right-justified in a field of width 32.
    """
    s = f"{x:.16E}"          # e.g. '1.2345678901234567E+00'
    mantissa, exp = s.split('E')
    sign = exp[0]
    digits = exp[1:].lstrip('0') or '0'
    exp3 = f"{sign}{int(digits):03d}"
    out = f"{mantissa}E{exp3}"
    return out.rjust(32)


def write_track(x: Track, state: EEPState) -> None:
    """Write ``x`` to ``x.filename`` in the Fortran ``.eep`` text format."""
    ncol = x.ncol
    if x.write_distance:
        ncol += 1
    if x.has_phase:
        have_phase = 'YES'.rjust(8)
        ncol += 1
    else:
        have_phase = ' NO'.rjust(8)

    dashes = '# ' + '-' * 86  # 88-char rule
    lines: List[str] = []

    lines.append(f"# MIST version number  = {x.version_string:<8.8}")
    lines.append(f"# MESA revision number = {x.MESA_revision_number:8d}")
    lines.append(dashes)
    lines.append('#  Yinit        Zinit   [Fe/H]   [a/Fe]  v/vcrit'.ljust(88))
    # '# ', f6.4, 1p1e13.5, 0p3f9.2
    yinit = f"{x.initial_Y:6.4f}"
    zinit = _fmt_1p_e13_5(x.initial_Z)
    rest = f"{x.Fe_div_H:9.2f}{x.alpha_div_Fe:9.2f}{x.v_div_vcrit:9.2f}"
    lines.append(f"# {yinit}{zinit}{rest}")
    lines.append(dashes)

    # column header describing the following two summary values
    lines.append(f"#{' ':1}{'initial_mass':>16}"
                 f"{'N_pts':>8}{'N_EEP':>8}{'N_col':>8}{'phase':>8}"
                 f"{' ':2}{'type':>10}")
    imass = _fmt_1p_e16_10(x.initial_mass)
    lines.append(f"#{' ':1}{imass}"
                 f"{x.ntrack:8d}{x.neep:8d}{ncol:8d}{have_phase}"
                 f"{' ':2}{STAR_LABEL[x.star_type - 1]:>10}")

    # EEPs line: '# EEPs: ' then the 1-based EEP locations, i8 each
    eeps = ''.join(f"{int(e):8d}" for e in x.eep)
    lines.append(f"{'# EEPs: ':<8}{eeps}")
    lines.append(dashes)

    # column-number row: each number right-justified in width 32 (27x + i5)
    col_numbers = '#' + ''.join(f"{' ':27}{k:5d}" for k in range(1, ncol + 1))
    lines.append(col_numbers)

    # column-name row
    names = [c.name for c in x.cols]
    if x.has_phase and x.write_distance:
        names = names + ['phase', 'distance']
    elif x.has_phase:
        names = names + ['phase']
    elif x.write_distance:
        names = names + ['distance']
    name_row = '#' + ''.join(f"{nm:>32}" for nm in names)
    lines.append(name_row)

    # data rows from eep(1) .. ntrack (1-based inclusive)
    start = int(x.eep[0])  # 1-based
    for j in range(start, x.ntrack + 1):  # inclusive
        col = j - 1
        vals = [x.tr[i, col] for i in range(x.ncol)]
        if x.has_phase:
            vals.append(float(x.phase[col]))
        if x.write_distance:
            vals.append(float(x.dist[col]))
        row = ' ' + ''.join(_fmt_es32_16e3(v) for v in vals)
        lines.append(row)

    with open(x.filename, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


def _fmt_1p_e13_5(x: float) -> str:
    """Fortran ``1p1e13.5`` -> one digit before point, 5 after, width 13."""
    s = f"{x:.5E}"                     # '1.42857E-02'
    mantissa, exp = s.split('E')
    sign = exp[0]
    val = int(exp[1:])
    return f"{mantissa}E{sign}{val:02d}".rjust(13)


def _fmt_1p_e16_10(x: float) -> str:
    """Fortran ``1p1e16.10`` -> one digit before point, 10 after, width 16."""
    s = f"{x:.10E}"
    mantissa, exp = s.split('E')
    sign = exp[0]
    val = int(exp[1:])
    return f"{mantissa}E{sign}{val:02d}".rjust(16)
