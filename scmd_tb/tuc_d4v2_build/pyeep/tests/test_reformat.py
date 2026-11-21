"""Tests for the MCDSEP -> MESA reformatter (``pyeep.reformat``).

Builds a synthetic MCDSEP file (one value per line), converts it, and checks
the result is valid fixed-width MESA that ``make_eep`` can read back.

Run directly (``python tests/test_reformat.py``) or under pytest.
"""

from __future__ import annotations

import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PY_ROOT = os.path.dirname(HERE)
sys.path.insert(0, PY_ROOT)
sys.path.insert(0, HERE)

from pyeep.reformat import (  # noqa: E402
    COLUMN_WIDTH,
    DEFAULT_COLUMNS,
    is_mesa_format,
    reformat_mesa_file,
    sniff_columns,
)

COLUMNS = list(DEFAULT_COLUMNS)


def write_mcdsep(path, nmodel, columns=COLUMNS, trailing_partial=0):
    """Write a synthetic MCDSEP file: 6 header lines then one value per line."""
    lines = [
        '  '.join(str(i) for i in range(1, 3)),
        '  '.join(['version_number', 'compiler']),
        '  '.join(['15140', 'gfortran']),
        '',
        '  '.join(str(i) for i in range(1, len(columns) + 1)),
        '  '.join(columns),
    ]
    for m in range(nmodel):
        for c in range(len(columns)):
            lines.append(f"   {float(m * 100 + c):.16E}")
    for c in range(trailing_partial):
        lines.append(f"   {float(c):.16E}")
    with open(path, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')
    return path


def read_rows(path):
    with open(path) as fh:
        lines = [ln.rstrip('\n') for ln in fh]
    return lines


def test_roundtrip():
    with tempfile.TemporaryDirectory() as tmp:
        src = write_mcdsep(os.path.join(tmp, 'raw'), nmodel=7)
        dst = os.path.join(tmp, 'out')
        n = reformat_mesa_file(src, dst)
        assert n == 7, n

        lines = read_rows(dst)
        # 6 header lines + 7 data rows
        assert len(lines) == 6 + 7, len(lines)

        # header layout matches what make_eep / input.format expect
        assert lines[1].split() == [
            'version_number', 'compiler', 'build', 'MESA_SDK_version',
            'math_backend', 'date', 'burn_min1', 'burn_min2']
        assert lines[3].strip() == ''
        assert lines[5].split() == COLUMNS

        # every data row has one field per column, and the values round-trip
        for j, row in enumerate(lines[6:]):
            fields = row.split()
            assert len(fields) == len(COLUMNS), (j, len(fields))
            vals = [float(x) for x in fields]
            assert vals == [float(j * 100 + c) for c in range(len(COLUMNS))]

        # fixed-width layout: width 40, single-space separator
        for row in lines[4:]:
            if not row.strip():
                continue
            assert len(row) == len(COLUMNS) * (COLUMN_WIDTH + 1) - 1, len(row)
        print("roundtrip OK")


def test_trailing_partial_model_dropped():
    """A run interrupted mid-model must not produce a short final row."""
    with tempfile.TemporaryDirectory() as tmp:
        src = write_mcdsep(os.path.join(tmp, 'raw'), nmodel=4, trailing_partial=5)
        dst = os.path.join(tmp, 'out')
        n = reformat_mesa_file(src, dst)
        assert n == 4, n
        lines = read_rows(dst)
        assert len(lines) == 6 + 4
        for row in lines[6:]:
            assert len(row.split()) == len(COLUMNS)
        print("trailing-partial OK")


def test_sniff_and_detect():
    with tempfile.TemporaryDirectory() as tmp:
        src = write_mcdsep(os.path.join(tmp, 'raw'), nmodel=3)
        with open(src) as fh:
            lines = [ln.rstrip('\n') for ln in fh]
        assert sniff_columns(lines) == COLUMNS

        # MCDSEP input is not MESA format; its output is
        assert not is_mesa_format(src)
        dst = os.path.join(tmp, 'out')
        reformat_mesa_file(src, dst)
        assert is_mesa_format(dst)
        print("sniff/detect OK")


def test_custom_columns():
    cols = ['star_age', 'log_L', 'log_Teff']
    with tempfile.TemporaryDirectory() as tmp:
        src = write_mcdsep(os.path.join(tmp, 'raw'), nmodel=5, columns=cols)
        dst = os.path.join(tmp, 'out')
        assert reformat_mesa_file(src, dst) == 5
        assert read_rows(dst)[5].split() == cols
        # explicit override wins over the file's own header
        dst2 = os.path.join(tmp, 'out2')
        assert reformat_mesa_file(src, dst2, columns=['a', 'b', 'c']) == 5
        assert read_rows(dst2)[5].split() == ['a', 'b', 'c']
        print("custom columns OK")


def test_no_complete_model_raises():
    with tempfile.TemporaryDirectory() as tmp:
        src = write_mcdsep(os.path.join(tmp, 'raw'), nmodel=0, trailing_partial=3)
        try:
            reformat_mesa_file(src, os.path.join(tmp, 'out'))
        except ValueError:
            print("empty-input raises OK")
        else:
            raise AssertionError("expected ValueError for a file with no models")


def test_make_eep_reads_reformatted():
    """The converted file must be readable by the make_eep history reader."""
    from pyeep.support import EEPState, Track, setup_columns, read_history_file

    with tempfile.TemporaryDirectory() as tmp:
        src = write_mcdsep(os.path.join(tmp, 'raw'), nmodel=40)
        dst = os.path.join(tmp, 'track')
        reformat_mesa_file(src, dst)

        cols_list = os.path.join(tmp, 'cols.list')
        with open(cols_list, 'w') as fh:
            fh.write('\n'.join('      ' + c for c in COLUMNS) + '\n')

        state = EEPState()
        setup_columns(cols_list, state)
        t = Track()
        t.filename = 'track'
        t.ncol = state.ncol
        t.neep = state.primary
        read_history_file(t, state, tmp)

        assert t.ntrack == 40, t.ntrack
        assert t.MESA_revision_number == 15140, t.MESA_revision_number
        assert t.tr.shape == (state.ncol, 40), t.tr.shape
        print("make_eep reads reformatted output OK")


if __name__ == '__main__':
    test_roundtrip()
    test_trailing_partial_model_dropped()
    test_sniff_and_detect()
    test_custom_columns()
    test_no_complete_model_raises()
    test_make_eep_reads_reformatted()
    print("\nall reformat tests passed")
