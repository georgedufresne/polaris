#!/usr/bin/env python3
"""Run the EEP pipeline on real MESA history files.

This is a convenience wrapper around :mod:`pyeep` for processing actual MESA
``history.data`` files without having to hand-write the two-part
``input.example`` file that :mod:`make_eep` expects.  It:

* accepts one or more history files (or ``LOGS`` directories containing a
  ``history.data``) on the command line,
* derives the ``history_columns.list`` directly from each file's own
  data-column header row (so whatever columns MESA wrote are used),
* reads ``&eep_controls`` from ``input.nml`` and the EEP intervals from
  ``input.eep`` if those files are present (falling back to sensible defaults),
* finds primary and secondary EEPs, assigns FSPS phases, and
* writes ``<name>.eep`` into an output directory (``python/output`` by default).

Examples
--------
Process a single history file::

    ./run_history.py /path/to/LOGS/history.data

Process every LOGS directory under a grid, naming outputs by mass::

    ./run_history.py ~/MESA_runs/*/LOGS -o output

Override composition metadata written into the .eep header::

    ./run_history.py history.data --initial-Y 0.27 --initial-Z 0.014
"""

from __future__ import annotations

import argparse
import copy
import os
import sys
from typing import List, Optional

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
from make_eep import read_namelist_group, apply_eep_controls  # noqa: E402


def resolve_history_path(arg: str) -> str:
    """Return a concrete history-file path from a file or LOGS directory arg."""
    if os.path.isdir(arg):
        candidate = os.path.join(arg, 'history.data')
        if os.path.exists(candidate):
            return candidate
        raise FileNotFoundError(f"no history.data found in directory {arg}")
    if os.path.exists(arg):
        return arg
    raise FileNotFoundError(f"history file not found: {arg}")


def derive_columns_list(history_path: str, dest_dir: str, name: str) -> str:
    """Write a history_columns.list from the file's own data-column header.

    MESA history layout puts the data-column *names* on line 6 (1-based).  We
    copy those verbatim so every column present in the file is carried through
    to the ``.eep`` output and is available for EEP identification.

    The file name is qualified with ``name`` (the track's output base name)
    because columns may differ between runs: a single shared
    ``history_columns.list`` would be a race when tracks are processed in
    parallel, and one worker could read another's column set.
    """
    with open(history_path) as fh:
        lines = fh.read().splitlines()
    data_col_names = lines[5].split()
    os.makedirs(dest_dir, exist_ok=True)
    dest = os.path.join(dest_dir, f'history_columns.{name}.list')
    with open(dest, 'w') as out:
        out.write("! columns derived from the MESA history file\n")
        out.write(f"! source: {history_path}\n\n")
        for name_ in data_col_names:
            out.write(f"      {name_}\n")
    return dest


def output_name(history_path: str) -> str:
    """Pick a sensible .eep base name for a given history file.

    Uses the parent directory name when the file is the generic ``history.data``
    (e.g. ``.../00080M/LOGS/history.data`` -> ``00080M``); otherwise uses the
    file's own stem.
    """
    base = os.path.basename(history_path)
    stem, _ = os.path.splitext(base)
    if stem == 'history':
        parent = os.path.dirname(os.path.abspath(history_path))
        name = os.path.basename(parent)
        if name.lower() == 'logs':
            name = os.path.basename(os.path.dirname(parent))
        if name:
            return name
    return stem


def build_state(args: argparse.Namespace) -> EEPState:
    """Create the shared EEPState, honouring input.nml / input.eep if present."""
    state = EEPState()

    nml = read_namelist_group(args.nml, 'eep_controls')
    flags = apply_eep_controls(nml, state)
    # allow an explicit --eep-input to win over the namelist value
    eep_input_file = args.eep_input or flags['eep_input_file']

    try:
        set_eep_interval(eep_input_file, state)
    except (OSError, ValueError) as exc:
        print(f"  note: could not read EEP intervals from {eep_input_file} "
              f"({exc}); using defaults", file=sys.stderr)
        set_default_eep_interval(state)

    state._flags = flags  # stash driver flags for run loop
    return state


def process_one(history_path: str, state: EEPState, args: argparse.Namespace) -> Optional[str]:
    """Run the full EEP pipeline on one history file; return the output path."""
    name = output_name(history_path)
    history_dir = os.path.dirname(os.path.abspath(history_path))
    history_file = os.path.basename(history_path)

    # per-file column list (columns can differ between runs)
    cols_list = derive_columns_list(history_path, args.output, name)
    setup_columns(cols_list, state)

    t = Track()
    t.filename = history_file
    t.ncol = state.ncol
    t.neep = state.primary

    read_history_file(t, state, history_dir)

    flags = state._flags
    t.he_star = flags['do_he_star']
    t.initial_Y = args.initial_Y
    t.initial_Z = args.initial_Z
    t.Fe_div_H = args.Fe_div_H
    t.alpha_div_Fe = args.alpha_div_Fe
    t.v_div_vcrit = args.v_div_vcrit
    t.version_string = args.version
    t.merger = False

    primary_eep(t, state)
    eep_str = ''.join(f"{int(e):8d}" for e in t.eep)
    print(f"{name}: initial_mass={t.initial_mass:.4f}  "
          f"ntrack={t.ntrack}  MESA={t.MESA_revision_number}")
    print(f"  primary EEPs:{eep_str}")

    if bool((t.eep == 0).all()):
        print(f"  PROBLEM: no EEPs defined for {name}; skipping")
        return None

    s = secondary_eep(t, state)
    if flags['do_phases']:
        set_track_phase(s, state)

    os.makedirs(args.output, exist_ok=True)
    out_path = os.path.join(args.output, f"{name}.eep")
    s.filename = out_path
    write_track(s, state)
    print(f"  wrote {out_path}  (N_pts={s.ntrack}, N_EEP={s.neep})")
    return out_path


class _HistoryJob:
    """One history file's work, in a picklable form for a worker pool.

    Each call gets its own copy of ``state`` because ``process_one`` mutates it
    via ``setup_columns`` (columns may differ between runs).  Copying keeps the
    job safe under the thread backend as well as the process backend.
    """

    def __init__(self, state: EEPState, args: argparse.Namespace):
        self.state = state
        self.args = args

    def __call__(self, history_path: str) -> Optional[str]:
        return process_one(history_path, copy.deepcopy(self.state), self.args)


def main(argv: Optional[List[str]] = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    parser = argparse.ArgumentParser(
        description="Run the pyeep EEP pipeline on real MESA history files.")
    parser.add_argument('history', nargs='+',
                        help="MESA history.data file(s) or LOGS directory(ies)")
    parser.add_argument('-o', '--output', default='output',
                        help="output directory for .eep files (default: output)")
    parser.add_argument('--nml', default='input.nml',
                        help="namelist file with &eep_controls (default: input.nml)")
    parser.add_argument('--eep-input', default=None,
                        help="EEP-interval file (overrides eep_input_file in the namelist)")
    parser.add_argument('--version', default='pyeep',
                        help="version string written to the .eep header (max 8 chars)")
    parser.add_argument('--initial-Y', type=float, default=0.0, dest='initial_Y')
    parser.add_argument('--initial-Z', type=float, default=0.0, dest='initial_Z')
    parser.add_argument('--Fe-div-H', type=float, default=0.0, dest='Fe_div_H')
    parser.add_argument('--alpha-div-Fe', type=float, default=0.0, dest='alpha_div_Fe')
    parser.add_argument('--v-div-vcrit', type=float, default=0.0, dest='v_div_vcrit')
    add_parallel_args(parser, default=1)
    args = parser.parse_args(argv)
    args.version = args.version[:8]

    state = build_state(args)

    paths: List[str] = []
    failed = 0
    for entry in args.history:
        try:
            paths.append(resolve_history_path(entry))
        except FileNotFoundError as exc:
            print(f"skip: {exc}", file=sys.stderr)
            failed += 1

    if not paths:
        print(f"\nDone: 0 written, {failed} skipped/failed", file=sys.stderr)
        return 1

    os.makedirs(args.output, exist_ok=True)

    workers = resolve_workers(args.jobs, len(paths))
    if workers > 1:
        print(f"processing {len(paths)} tracks on {workers} {args.backend}",
              file=sys.stderr)

    written = 0
    job = _HistoryJob(state, args)
    for path, result in map_ordered(job, paths, workers, args.backend):
        if isinstance(result, Exception):
            print(f"error processing {path}: {result}", file=sys.stderr)
            failed += 1
        elif result:
            written += 1
        else:
            failed += 1

    print(f"\nDone: {written} written, {failed} skipped/failed "
          f"-> {os.path.abspath(args.output)}")
    return 0 if written > 0 else 1


if __name__ == '__main__':
    raise SystemExit(main())
