"""End-to-end test for the Python make_eep port.

Builds a synthetic MESA history file plus the required input files in a
temporary directory, runs the driver, and verifies the resulting ``.eep`` file
is well-formed (correct header fields and numeric columns).

Run directly (``python tests/test_make_eep.py``) or under pytest.
"""

from __future__ import annotations

import os
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
PY_ROOT = os.path.dirname(HERE)
sys.path.insert(0, PY_ROOT)
sys.path.insert(0, HERE)

import make_eep  # noqa: E402
from make_synthetic_track import (  # noqa: E402
    write_history_file,
    write_history_columns_list,
)

INPUT_EEP = """#format i3, i5
9
  1  200
  2  150
  3  100
  4  150
  5   25
  6   75
  7  100
  8  600
  9  300
"""

INPUT_NML = """&eep_controls
   eep_input_file = 'input.eep'
   do_phases = .true.
   do_he_star = .false.
   make_bin_tracks = .false.
   center_gamma_limit = 100
   center_carbon_limit = 0.05
   log_center_T_limit = 8.5
   high_mass_limit = 10.0
   very_low_mass_limit = 0.7
   weight_center_rho_T_by_Xc = .false.
   Tc_scale = 1d-2
   Rhoc_scale = 1d-2
   age_scale = 1d0
   Teff_scale = 2.0d0
   logL_scale = 0.125d0
/
"""

INPUT_FORMAT = """# header on line 1, data on line 2
40
40
1
"""


def _build_case(root: str) -> str:
    tracks_dir = os.path.join(root, 'tracks')
    eeps_dir = os.path.join(root, 'eeps')
    os.makedirs(tracks_dir, exist_ok=True)
    os.makedirs(eeps_dir, exist_ok=True)

    write_history_file(os.path.join(tracks_dir, '00100M.track'),
                       initial_mass=1.0)
    cols_list = os.path.join(root, 'history_columns.list')
    write_history_columns_list(cols_list)

    with open(os.path.join(root, 'input.eep'), 'w') as fh:
        fh.write(INPUT_EEP)
    with open(os.path.join(root, 'input.nml'), 'w') as fh:
        fh.write(INPUT_NML)
    with open(os.path.join(root, 'input.format'), 'w') as fh:
        fh.write(INPUT_FORMAT)

    input_example = os.path.join(root, 'input.example')
    with open(input_example, 'w') as fh:
        fh.write("#version string, max 8 characters\n")
        fh.write("testver\n")
        fh.write("#initial Y, Z, [Fe/H], [alpha/Fe], v/vcrit\n")
        fh.write("   0.2703 1.42857e-02   0.00        0.00     0.4\n")
        fh.write("#data directories\n")
        fh.write(f"{tracks_dir}\n")
        fh.write(f"{eeps_dir}\n")
        fh.write(f"{os.path.join(root, 'iso')}\n")
        fh.write("#read history_columns\n")
        fh.write(f"{cols_list}\n")
        fh.write("#specify tracks\n")
        fh.write("1\n")
        fh.write("00100M.track\n")
    return input_example, eeps_dir


def _parse_eep(path: str):
    with open(path) as fh:
        lines = fh.read().splitlines()
    header = {}
    header['version'] = lines[0].split('=', 1)[1].strip()
    header['mesa_rev'] = int(lines[1].split('=', 1)[1])
    # summary values line is the 8th line (index 7)
    summary = lines[7].split()
    # summary: '#', initial_mass, N_pts, N_EEP, N_col, phase(YES/NO), type...
    header['initial_mass'] = float(summary[1])
    header['n_pts'] = int(summary[2])
    header['n_eep'] = int(summary[3])
    header['n_col'] = int(summary[4])
    header['phase'] = summary[5]
    # EEPs line (index 8) begins with '# EEPs:'
    eeps = [int(x) for x in lines[8].replace('# EEPs:', '').split()]
    header['eeps'] = eeps
    # find first data line (does not start with '#')
    data = [ln for ln in lines if ln and not ln.lstrip().startswith('#')]
    return header, data


def test_make_eep_end_to_end():
    with tempfile.TemporaryDirectory() as root:
        cwd = os.getcwd()
        input_example, eeps_dir = _build_case(root)
        try:
            os.chdir(root)
            rc = make_eep.run('input.example')
        finally:
            os.chdir(cwd)

        assert rc == 0
        out = os.path.join(eeps_dir, '00100M.track.eep')
        assert os.path.exists(out), "eep output file was not written"

        header, data = _parse_eep(out)

        assert header['version'] == 'testver'
        assert header['mesa_rev'] == 12000
        assert abs(header['initial_mass'] - 1.0) < 1e-6
        assert header['phase'] == 'YES'

        # at least the H-burning EEPs (PreMS, ZAMS, TAMS, TAMS) must be found
        nonzero = [e for e in header['eeps'] if e > 0]
        assert len(nonzero) >= 4, f"too few EEPs found: {header['eeps']}"

        # first EEP must be at 1 (PreMS_fudge)
        assert header['eeps'][0] == 1

        # number of data rows equals N_pts
        assert len(data) == header['n_pts'], (
            f"data rows {len(data)} != N_pts {header['n_pts']}")

        # each data row has n_col numeric fields (history cols + phase)
        ncols = header['n_col']
        for row in data[:5] + data[-5:]:
            fields = row.split()
            assert len(fields) == ncols, (
                f"row has {len(fields)} fields, expected {ncols}")
            # all parse as floats
            [float(x) for x in fields]

        print("EEPs found:", header['eeps'])
        print(f"N_pts={header['n_pts']} N_EEP={header['n_eep']} "
              f"N_col={header['n_col']}")
        print("OK")


def _build_real_case(root: str, history_data: str):
    """Build an input case around a *real* MESA history.data file."""
    tracks_dir = os.path.join(root, 'tracks')
    eeps_dir = os.path.join(root, 'eeps')
    os.makedirs(tracks_dir, exist_ok=True)
    os.makedirs(eeps_dir, exist_ok=True)

    import shutil
    shutil.copy(history_data, os.path.join(tracks_dir, '00080M.track'))

    # derive history_columns.list from the file's own data-column header (line 6)
    with open(history_data) as fh:
        lines = fh.read().splitlines()
    data_col_names = lines[5].split()
    cols_list = os.path.join(root, 'history_columns.list')
    with open(cols_list, 'w') as fh:
        fh.write("! columns derived from the real MESA history file\n\n")
        for c in data_col_names:
            fh.write(f"      {c}\n")

    with open(os.path.join(root, 'input.eep'), 'w') as fh:
        fh.write(INPUT_EEP)
    with open(os.path.join(root, 'input.nml'), 'w') as fh:
        fh.write(INPUT_NML)
    with open(os.path.join(root, 'input.format'), 'w') as fh:
        fh.write(INPUT_FORMAT)

    input_example = os.path.join(root, 'input.example')
    with open(input_example, 'w') as fh:
        fh.write("#version string\n")
        fh.write("mist08\n")
        fh.write("#initial Y, Z, [Fe/H], [alpha/Fe], v/vcrit\n")
        fh.write("   0.2703 1.42857e-02   0.00        0.00     0.4\n")
        fh.write("#data directories\n")
        fh.write(f"{tracks_dir}\n{eeps_dir}\n{os.path.join(root, 'iso')}\n")
        fh.write("#read history_columns\n")
        fh.write(f"{cols_list}\n")
        fh.write("#specify tracks\n1\n00080M.track\n")
    return input_example, eeps_dir


def test_make_eep_real_data():
    """Optional regression test against a real MESA history file.

    Enabled by setting the ``PYEEP_REAL_HISTORY`` environment variable to the
    path of a MESA ``history.data`` file; otherwise it is skipped so the suite
    remains self-contained.
    """
    history_data = os.environ.get('PYEEP_REAL_HISTORY')
    if not history_data or not os.path.exists(history_data):
        print("SKIP real-data test (set PYEEP_REAL_HISTORY to enable)")
        return

    with tempfile.TemporaryDirectory() as root:
        cwd = os.getcwd()
        input_example, eeps_dir = _build_real_case(root, history_data)
        try:
            os.chdir(root)
            rc = make_eep.run('input.example')
        finally:
            os.chdir(cwd)

        assert rc == 0
        out = os.path.join(eeps_dir, '00080M.track.eep')
        assert os.path.exists(out), "eep output not written for real data"

        header, data = _parse_eep(out)
        # first EEP is always PreMS at index 1
        assert header['eeps'][0] == 1
        # every data row parses and has N_col fields
        ncols = header['n_col']
        for row in data[:3] + data[-3:]:
            fields = row.split()
            assert len(fields) == ncols
            [float(x) for x in fields]
        assert len(data) == header['n_pts']
        print("REAL EEPs found:", header['eeps'])
        print(f"REAL N_pts={header['n_pts']} N_EEP={header['n_eep']} "
              f"N_col={header['n_col']}  OK")


if __name__ == '__main__':
    test_make_eep_end_to_end()
    test_make_eep_real_data()
