"""Convert MCDSEP-style track files into the fixed-width MESA history format.

Some track generators (e.g. the Dartmouth/MCDSEP pipeline) write history data
with **one value per line** rather than one model per line: after a six-line
header the file is a flat stream of numbers, 22 of them (one per column) for
each model.  ``make_eep`` — both the Fortran original and this port — expects
standard MESA layout, where each model occupies a single row.

This module rewrites the former into the latter.  The output uses a fixed
column width of 40 characters separated by a single space, which matches the
``input.format`` defaults shipped with ``iso``::

    40      # head  -- header field width
    40      # main  -- data field width
    1       # xtra  -- separator between fields

Use it from the command line via ``reformat_mesa.py``, or programmatically::

    from pyeep.reformat import reformat_mesa_file
    reformat_mesa_file('m0700.feh031.20001', 'm0700.feh031.20001_mesa23')
"""

from __future__ import annotations

import os
from typing import List, Optional, Sequence

#: Field width for every column, header and data alike.
COLUMN_WIDTH = 40

#: Number of header lines in an MCDSEP file before the flat value stream.
RAW_HEADER_LINES = 6

#: Column layout assumed when a file carries no usable column-name header.
DEFAULT_COLUMNS: Sequence[str] = (
    'star_age', 'star_mass', 'star_mdot', 'he_core_mass', 'c_core_mass',
    'log_L', 'log_LH', 'log_LHe', 'log_Teff', 'log_R', 'log_g',
    'surface_h1', 'surface_he3', 'surface_he4', 'surface_c12', 'surface_o16',
    'log_center_T', 'log_center_Rho', 'center_gamma', 'center_h1',
    'center_he4', 'center_c12',
)

#: Names of the eight fields written into the MESA file header.
_HEADER_NAMES: Sequence[str] = (
    'version_number', 'compiler', 'build', 'MESA_SDK_version',
    'math_backend', 'date', 'burn_min1', 'burn_min2',
)

#: Default values for those fields, used when the source file has none.
_HEADER_VALUES: Sequence[str] = (
    '15140', '"gfortran"', '"10.2.0"', '"x86_64-linux-20.12.1"',
    '"CRMATH"', '"20210328"', '5.0000000000000000E+001',
    '1.0000000000000000E+003',
)


def _pad(field: str) -> str:
    """Right-justify ``field`` in a fixed-width column, truncating if needed."""
    return f"{field:>{COLUMN_WIDTH}}"[:COLUMN_WIDTH]


def _row(fields: Sequence[str]) -> str:
    return ' '.join(_pad(f) for f in fields)


def _read_lines(path: str, encoding: Optional[str] = None) -> List[str]:
    """Read ``path`` as text, tolerating unknown / non-UTF-8 encodings.

    When ``encoding`` is not given, ``chardet`` is used if installed, then
    UTF-8, then latin-1 (which cannot fail).  Tabs are normalised to spaces.
    """
    encodings: List[str] = []
    if encoding:
        encodings.append(encoding)
    else:
        try:
            import chardet  # optional dependency
            with open(path, 'rb') as fh:
                guess = chardet.detect(fh.read(10000)).get('encoding')
            if guess:
                encodings.append(guess)
        except ImportError:
            pass
        encodings.append('utf-8')
    encodings.append('latin-1')  # single-byte fallback, never raises

    last: Optional[Exception] = None
    for enc in encodings:
        try:
            with open(path, 'r', encoding=enc) as fh:
                return [ln.replace('\t', ' ').rstrip('\n') for ln in fh]
        except (UnicodeDecodeError, LookupError) as exc:
            last = exc
    raise UnicodeDecodeError(  # pragma: no cover - latin-1 makes this unreachable
        'pyeep', b'', 0, 1, f'could not decode {path}: {last}')


def sniff_columns(lines: Sequence[str]) -> Optional[List[str]]:
    """Return the column names from an MCDSEP header, or ``None`` if absent.

    The sixth line of an MCDSEP file holds the whitespace-separated column
    names.  It is only accepted when it looks like a name row (at least two
    entries, none of them parsable as a number).
    """
    if len(lines) < RAW_HEADER_LINES:
        return None
    names = lines[RAW_HEADER_LINES - 1].split()
    if len(names) < 2:
        return None
    for name in names:
        try:
            float(name)
        except ValueError:
            continue
        return None  # numeric entry -> this is data, not a name row
    return names


def reformat_mesa_file(input_file: str,
                       output_file: str,
                       columns: Optional[Sequence[str]] = None,
                       encoding: Optional[str] = None,
                       header_values: Optional[Sequence[str]] = None) -> int:
    """Rewrite ``input_file`` as a MESA-format history file at ``output_file``.

    Parameters
    ----------
    input_file, output_file
        Source MCDSEP file and destination path.
    columns
        Column names.  Defaults to the names found in the source header,
        falling back to :data:`DEFAULT_COLUMNS`.
    encoding
        Force a text encoding instead of detecting one.
    header_values
        The eight MESA header values.  Defaults to :data:`_HEADER_VALUES`.

    Returns
    -------
    int
        The number of complete models (data rows) written.

    Raises
    ------
    ValueError
        If the file holds no complete model.
    """
    lines = _read_lines(input_file, encoding)

    if columns is None:
        columns = sniff_columns(lines) or list(DEFAULT_COLUMNS)
    ncol = len(columns)
    if ncol < 1:
        raise ValueError(f"{input_file}: no columns to write")

    values = header_values if header_values is not None else _HEADER_VALUES
    if len(values) != len(_HEADER_NAMES):
        raise ValueError(f"header_values must have {len(_HEADER_NAMES)} entries")

    data = lines[RAW_HEADER_LINES:]

    out: List[str] = [
        _row([str(i) for i in range(1, len(_HEADER_NAMES) + 1)]),
        _row(_HEADER_NAMES),
        _row(values),
        '',
        _row([str(i) for i in range(1, ncol + 1)]),
        _row(columns),
    ]

    nrow = 0
    for i in range(0, len(data), ncol):
        chunk = data[i:i + ncol]
        if len(chunk) < ncol:
            break  # trailing partial model -- a truncated/interrupted run
        out.append(_row([f.strip() for f in chunk]))
        nrow += 1

    if nrow == 0:
        raise ValueError(
            f"{input_file}: no complete models "
            f"({len(data)} data lines, need a multiple of {ncol})")

    parent = os.path.dirname(output_file)
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(output_file, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(out) + '\n')

    return nrow


def is_mesa_format(path: str, encoding: Optional[str] = None) -> bool:
    """Return True if ``path`` already looks like a row-per-model MESA file.

    Used to skip files that need no conversion: in MESA layout the first data
    line carries one value per column, whereas MCDSEP files carry exactly one.
    """
    try:
        lines = _read_lines(path, encoding)
    except OSError:
        return False
    if len(lines) <= RAW_HEADER_LINES:
        return False
    first = lines[RAW_HEADER_LINES].split()
    return len(first) > 1
