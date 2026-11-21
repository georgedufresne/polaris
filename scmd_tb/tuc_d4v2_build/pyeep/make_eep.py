#!/usr/bin/env python3
"""make_eep: convert MESA history files into EEP-based tracks.

Python port of ``src/make_eep.f90``.  Usage::

    ./make_eep.py <input>

where ``<input>`` is the two-part input file (e.g. ``input.example``).  As in
the Fortran version, three auxiliary files are read from the current working
directory:

* ``input.nml``    -- the ``&eep_controls`` namelist (tunable parameters),
* ``input.format`` -- the history-file column widths (kept for reference),
* the EEP-interval file named by ``eep_input_file`` (default ``input.eep``).

For each track listed in ``<input>`` the program reads the MESA history file,
identifies primary and secondary EEPs, optionally assigns FSPS phases and writes
an ``<name>.eep`` file into the eeps directory.
"""

from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional

# allow running directly from the python/ directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pyeep.support import (  # noqa: E402
    EEPState,
    Track,
    setup_columns,
    set_eep_interval,
    set_default_eep_interval,
    read_history_file,
    write_track,
)
from pyeep.eep import primary_eep, secondary_eep  # noqa: E402
from pyeep.phases import set_track_phase  # noqa: E402
from pyeep.parallel import add_parallel_args, map_ordered, resolve_workers  # noqa: E402


# ---------------------------------------------------------------------------
# minimal Fortran-namelist reader (only what eep_controls needs)
# ---------------------------------------------------------------------------

def _convert_value(raw: str):
    """Convert a Fortran namelist scalar literal to a Python value."""
    s = raw.strip().rstrip(',').strip()
    if len(s) >= 2 and s[0] == "'" and s[-1] == "'":
        return s[1:-1]
    if len(s) >= 2 and s[0] == '"' and s[-1] == '"':
        return s[1:-1]
    low = s.lower()
    if low in ('.true.', 't', '.t.'):
        return True
    if low in ('.false.', 'f', '.f.'):
        return False
    # numeric: Fortran uses d/D for double-precision exponents
    num = low.replace('d', 'e')
    try:
        if any(ch in num for ch in ('.', 'e')):
            return float(num)
        return int(num)
    except ValueError:
        return s


def read_namelist_group(path: str, group: str) -> Dict[str, object]:
    """Return a dict of key->value for the namelist ``&group`` in ``path``.

    Comment lines (starting with ``!``) and inline ``!`` comments are ignored.
    Only simple ``key = value`` scalar assignments are supported, which covers
    the ``eep_controls`` group.
    """
    result: Dict[str, object] = {}
    if not os.path.exists(path):
        return result
    in_group = False
    with open(path, 'r') as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith('!'):
                continue
            # strip trailing inline comment
            bang = line.find('!')
            if bang >= 0:
                line = line[:bang].strip()
            if not line:
                continue
            if not in_group:
                if line.lower().startswith('&' + group.lower()):
                    in_group = True
                continue
            if line.startswith('/'):
                break
            if line.startswith('&'):
                break
            if '=' in line:
                key, val = line.split('=', 1)
                result[key.strip().lower()] = _convert_value(val)
    return result


# ---------------------------------------------------------------------------
# input.example reader (read_input)
# ---------------------------------------------------------------------------

@dataclass
class InputSpec:
    version_string: str = 'unknown'
    initial_Y: float = 0.0
    initial_Z: float = 0.0
    Fe_div_H: float = 0.0
    alpha_div_Fe: float = 0.0
    v_div_vcrit: float = 0.0
    history_dir: str = ''
    eep_dir: str = ''
    iso_dir: str = ''
    history_columns_list: str = ''
    history_files: List[str] = field(default_factory=list)


def read_input_file(path: str) -> InputSpec:
    """Parse the track section of the two-part input file (input.example)."""
    with open(path, 'r') as fh:
        lines = [ln.rstrip('\n') for ln in fh]
    it = iter(lines)

    def nextline() -> str:
        return next(it)

    spec = InputSpec()
    nextline()                                   # comment
    spec.version_string = nextline()[:8]         # a8
    nextline()                                   # comment
    parts = nextline().split()
    (spec.initial_Y, spec.initial_Z, spec.Fe_div_H,
     spec.alpha_div_Fe, spec.v_div_vcrit) = (float(p) for p in parts[:5])
    nextline()                                   # comment
    spec.history_dir = nextline().strip()
    spec.eep_dir = nextline().strip()
    spec.iso_dir = nextline().strip()
    nextline()                                   # comment
    spec.history_columns_list = nextline().strip()
    nextline()                                   # comment
    num = int(nextline().split()[0])
    for _ in range(num):
        spec.history_files.append(nextline().strip())
    return spec


def read_format_file(path: str, state: EEPState) -> None:
    """Read input.format (head, main, xtra).  Kept for parity; tolerant of absence."""
    if not os.path.exists(path):
        return
    with open(path, 'r') as fh:
        lines = [ln for ln in fh]
    try:
        state.head = int(lines[1].split()[0])
        state.main = int(lines[2].split()[0])
        state.xtra = int(lines[3].split()[0])
    except (IndexError, ValueError):
        pass


# ---------------------------------------------------------------------------
# apply the eep_controls namelist to the shared state
# ---------------------------------------------------------------------------

def apply_eep_controls(nml: Dict[str, object], state: EEPState) -> dict:
    """Copy recognised eep_controls values into ``state``; return driver flags."""
    def get(name, default):
        return nml.get(name, default)

    state.center_gamma_limit = float(get('center_gamma_limit', state.center_gamma_limit))
    state.center_carbon_limit = float(get('center_carbon_limit', state.center_carbon_limit))
    state.log_center_T_limit = float(get('log_center_t_limit', state.log_center_T_limit))
    state.high_mass_limit = float(get('high_mass_limit', state.high_mass_limit))
    state.very_low_mass_limit = float(get('very_low_mass_limit', state.very_low_mass_limit))
    state.weight_center_rho_T_by_Xc = bool(get('weight_center_rho_t_by_xc', state.weight_center_rho_T_by_Xc))
    state.Teff_scale = float(get('teff_scale', state.Teff_scale))
    state.logL_scale = float(get('logl_scale', state.logL_scale))
    state.age_scale = float(get('age_scale', state.age_scale))
    state.Tc_scale = float(get('tc_scale', state.Tc_scale))
    state.Rhoc_scale = float(get('rhoc_scale', state.Rhoc_scale))
    state.IAMS_Xc = float(get('iams_xc', state.IAMS_Xc))
    state.RGBTip_require_He_ignition = bool(get('rgbtip_require_he_ignition',
                                                state.RGBTip_require_He_ignition))
    state.He_ignition_log_ratio = float(get('he_ignition_log_ratio', state.He_ignition_log_ratio))
    state.RGBTip_take_later = bool(get('rgbtip_take_later', state.RGBTip_take_later))
    state.add_RGBBase = bool(get('add_rgbbase', state.add_RGBBase))
    state.VLM_fallback_min_age = float(get('vlm_fallback_min_age', state.VLM_fallback_min_age))
    state.merger_fallback = bool(get('merger_fallback', state.merger_fallback))
    state.Xc_HeL_ladder = bool(get('xc_hel_ladder', state.Xc_HeL_ladder))

    return {
        'do_phases': bool(get('do_phases', True)),
        'do_he_star': bool(get('do_he_star', False)),
        'eep_input_file': str(get('eep_input_file', 'input.eep')),
    }


# ---------------------------------------------------------------------------
# main driver (make_eeps)
# ---------------------------------------------------------------------------

class _TrackJob:
    """One track's work, in a form that can be shipped to a worker process.

    Holds only picklable data (the state dataclass, the input spec fields and
    the filename), so the same callable works inline or in a pool.
    """

    def __init__(self, state: EEPState, spec: 'InputSpec', flags: dict):
        self.state = state
        self.spec = spec
        self.flags = flags

    def __call__(self, filename: str) -> bool:
        state, spec, flags = self.state, self.spec, self.flags

        t = Track()
        t.filename = filename
        t.ncol = state.ncol
        t.neep = state.primary

        try:
            read_history_file(t, state, spec.history_dir)
        except OSError as exc:
            print(f"make_eep: problem reading {filename}: {exc}", file=sys.stderr)
            return False

        t.he_star = flags['do_he_star']
        t.initial_Y = spec.initial_Y
        t.initial_Z = spec.initial_Z
        t.Fe_div_H = spec.Fe_div_H
        t.alpha_div_Fe = spec.alpha_div_Fe
        t.v_div_vcrit = spec.v_div_vcrit
        t.version_string = spec.version_string
        t.merger = False

        print(f"{t.filename} {t.neep} {t.MESA_revision_number}")

        primary_eep(t, state)
        print(''.join(f"{int(e):8d}" for e in t.eep))

        if bool((t.eep == 0).all()):
            print(" PROBLEM WITH TRACK: NO EEPS DEFINED ")
            return False

        s = secondary_eep(t, state)
        if flags['do_phases']:
            set_track_phase(s, state)
        s.filename = os.path.join(spec.eep_dir, f"{s.filename}.eep")
        os.makedirs(os.path.dirname(s.filename) or '.', exist_ok=True)
        write_track(s, state)
        return True


def run(input_file: str,
        nml_file: str = 'input.nml',
        format_file: str = 'input.format',
        jobs: Optional[int] = 1,
        backend: str = 'processes') -> int:
    state = EEPState()

    nml = read_namelist_group(nml_file, 'eep_controls')
    flags = apply_eep_controls(nml, state)

    spec = read_input_file(input_file)

    # EEP interval spec
    eep_input_file = flags['eep_input_file']
    try:
        set_eep_interval(eep_input_file, state)
    except (OSError, ValueError) as exc:
        print(f" make_eeps: problem reading eep input file ({exc})", file=sys.stderr)
        print("            setting default EEPs", file=sys.stderr)
        set_default_eep_interval(state)
    if state.add_RGBBase and state.primary < 10:
        print(f" make_eeps: add_RGBBase needs one more primary than {eep_input_file} "
              f"lists ({state.primary}); see run/input_dsee.eep", file=sys.stderr)
    if state.Xc_HeL_ladder:
        n_ladder = 3 + len(state.Xc_ladder) + len(state.HeL_ladder)
        if state.primary != n_ladder:
            print(f" make_eeps: xc_hel_ladder places {n_ladder} primary EEPs but "
                  f"{eep_input_file} lists {state.primary}; "
                  f"see run/input_dsee_ladder.eep", file=sys.stderr)
        if flags['do_phases']:
            print(" make_eeps: the ladder primaries are not the FSPS ones, so "
                  "do_phases is off", file=sys.stderr)
            flags['do_phases'] = False

    read_format_file(format_file, state)

    setup_columns(spec.history_columns_list, state)

    files = spec.history_files
    if not files:
        return 0

    # the eeps directory is shared by every worker; create it once up front so
    # the workers never race on it
    if spec.eep_dir:
        os.makedirs(spec.eep_dir, exist_ok=True)

    workers = resolve_workers(jobs, len(files))
    if workers > 1:
        print(f" make_eeps: processing {len(files)} tracks "
              f"on {workers} {backend}", file=sys.stderr)

    job = _TrackJob(state, spec, flags)
    for filename, result in map_ordered(job, files, workers, backend):
        if isinstance(result, Exception):
            print(f"make_eep: problem processing {filename}: {result}",
                  file=sys.stderr)

    return 0


def main(argv: Optional[List[str]] = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0] in ('-h', '--help'):
        print("   make_eeps                  ")
        print("   usage: ./make_eep.py <input> [-j N] [--backend B] ")
        if not argv:
            print("   no command line argument   ")
        print()
        print("   -j, --jobs N     parallel worker processes; 0 = auto-detect")
        print("                    (default: 1, i.e. serial)")
        print("   --backend B      'processes' (default) or 'threads'")
        return 0 if argv else 1

    parser = argparse.ArgumentParser(
        prog='make_eep.py', add_help=False,
        description="Convert MESA history files into EEP-based tracks.")
    parser.add_argument('input')
    add_parallel_args(parser, default=1)
    args = parser.parse_args(argv)

    return run(args.input, jobs=args.jobs, backend=args.backend)


if __name__ == '__main__':
    raise SystemExit(main())
