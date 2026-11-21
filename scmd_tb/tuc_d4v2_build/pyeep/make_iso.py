#!/usr/bin/env python3
"""make_iso: build isochrones from EEP-based tracks.

Python port of ``src/make_iso.f90``.  Usage::

    ./make_iso.py <input> [<input> ...] [-j N]

``<input>`` is the same two-part input file ``make_eep.py`` reads; this program
uses the *isochrone* half (output filename, age list, optional ``double``) and
reads the ``.eep`` files that make_eep wrote into the eeps directory.  As in the
Fortran, ``input.nml`` (the ``&iso_controls`` namelist) and the history column
list are read from the current working directory.

Several input files may be given -- one per grid -- and ``-j`` runs them on
that many processes.  Output per file is identical to a serial run.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import time
from dataclasses import dataclass, field
from typing import List, Optional

# allow running directly from any directory
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from make_eep import read_namelist_group  # noqa: E402
from pyeep.support import (  # noqa: E402
    AGE_SCALE_LINEAR,
    AGE_SCALE_LOG10,
    process_history_columns,
)
from pyeep.iso import (  # noqa: E402
    IsoControls,
    IsochroneGrid,
    MakeIsoError,
    read_eep,
    write_isochrones,
)
from pyeep.parallel import add_parallel_args, map_ordered, resolve_workers  # noqa: E402


# ---------------------------------------------------------------------------
# input file (read_iso_input)
# ---------------------------------------------------------------------------

@dataclass
class IsoInputSpec:
    version_string: str = ''
    initial_Y: float = 0.0
    initial_Z: float = 0.0
    Fe_div_H: float = 0.0
    alpha_div_Fe: float = 0.0
    v_div_vcrit: float = 0.0
    history_dir: str = ''
    eep_dir: str = ''
    iso_dir: str = ''
    history_columns_list: str = ''
    tracks: List[str] = field(default_factory=list)
    output: str = ''
    age_scale: int = AGE_SCALE_LINEAR
    ages: List[float] = field(default_factory=list)
    double_eep: bool = False


def _list_real(line: str) -> float:
    """First value of a list-directed real read (accepts Fortran d exponents)."""
    return float(re.split(r'[\s,]+', line.strip())[0].lower().replace('d', 'e'))


def read_iso_input(path: str) -> IsoInputSpec:
    """Parse both halves of the input file, positionally, as make_iso does."""
    with open(path, 'r') as fh:
        lines = [ln.rstrip('\n') for ln in fh]
    it = iter(lines)
    nxt = lambda: next(it)  # noqa: E731

    spec = IsoInputSpec()
    nxt()
    spec.version_string = nxt()[:8]
    nxt()
    parts = nxt().split()
    (spec.initial_Y, spec.initial_Z, spec.Fe_div_H,
     spec.alpha_div_Fe, spec.v_div_vcrit) = (float(p.lower().replace('d', 'e'))
                                            for p in parts[:5])
    nxt()
    spec.history_dir = nxt().strip()
    spec.eep_dir = nxt().strip()
    spec.iso_dir = nxt().strip()
    nxt()
    spec.history_columns_list = nxt().strip()
    nxt()
    ntrk = int(nxt().split()[0])
    spec.tracks = [nxt().strip() for _ in range(ntrk)]

    nxt()                                            # '#specify isochrones'
    # list-directed character read: first token, ended by blank, comma or /
    name = re.split(r'[\s,/]', nxt().strip().strip('\'"'))[0]
    spec.output = spec.iso_dir.rstrip() + '/' + name
    list_type = nxt()[:10].strip()
    age_type = nxt()[:10].strip()
    if age_type == 'linear':
        spec.age_scale = AGE_SCALE_LINEAR
    elif age_type == 'log10':
        spec.age_scale = AGE_SCALE_LOG10
    else:
        raise MakeIsoError(' make_iso: age scale must be given as "linear" or "log10"')
    niso = int(nxt().split()[0])
    if list_type == 'min_max':
        lo = _list_real(nxt())
        hi = _list_real(nxt())
        if hi < lo:
            raise MakeIsoError('  make_iso: max age < min age')
        step = (hi - lo) / float(niso - 1) if niso > 1 else 0.0
        spec.ages = [lo + step * float(i) for i in range(niso)]
    elif list_type == 'list':
        spec.ages = [_list_real(nxt()) for _ in range(niso)]
    else:
        raise MakeIsoError(' make_iso: ages must be given as "list" or "min_max"')

    try:
        spec.double_eep = nxt()[:6] == 'double'
    except StopIteration:
        spec.double_eep = False
    return spec


# ---------------------------------------------------------------------------
# reading the tracks
# ---------------------------------------------------------------------------

def load_tracks(spec: IsoInputSpec, debug: bool = False):
    """Read every listed .eep file, skipping missing ones, with make_iso's checks."""
    good = []
    for name in spec.tracks:
        path = f"{spec.eep_dir.rstrip()}/{name}.eep"
        try:
            t = read_eep(path)
        except OSError:
            print(f"  PROBLEM OPENING EEP FILE: {path}")
            continue
        if debug:
            print(f"{name:>50}{t.initial_mass:8.2f}"
                  + ''.join(f"{int(e):8d}" for e in t.eep))
        # the Fortran compares each track, from the third good one on,
        # against the first good track
        if len(good) > 1:
            ref = good[0]
            if t.initial_mass < ref.initial_mass:
                raise MakeIsoError(' make_iso: masses out of order')
            if t.has_phase != ref.has_phase:
                raise MakeIsoError(' make_iso: inconsistent phase info in tracks')
            if t.MESA_revision_number != ref.MESA_revision_number:
                raise MakeIsoError(' make_iso: inconsistent version number in tracks')
        good.append(t)
    return good


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

def make_isochrones(input_file: str, nml_file: str = 'input.nml',
                    write: bool = True):
    """Run make_iso for one input file; returns the list of isochrones."""
    nml = read_namelist_group(nml_file, 'iso_controls')
    if not os.path.exists(nml_file):
        raise MakeIsoError(f'make_iso: problem reading iso_controls namelist '
                           f'({nml_file} not found)')
    controls = IsoControls.from_namelist(nml)
    spec = read_iso_input(input_file)
    cols = process_history_columns(spec.history_columns_list)
    tracks = load_tracks(spec, controls.iso_debug)

    grid = IsochroneGrid(tracks, [c.name for c in cols], controls)
    grid.initial_Y, grid.initial_Z = spec.initial_Y, spec.initial_Z
    grid.Fe_div_H, grid.alpha_div_Fe = spec.Fe_div_H, spec.alpha_div_Fe
    grid.v_div_vcrit = spec.v_div_vcrit

    isos = [grid.isochrone(a, spec.age_scale, spec.double_eep)
            for a in spec.ages]
    if write:
        print(f"  isochrone output file = {spec.output}", file=sys.stderr)
        write_isochrones(spec.output, isos, spec.version_string,
                         tracks[0].MESA_revision_number, spec.initial_Y,
                         spec.initial_Z, spec.Fe_div_H, spec.alpha_div_Fe,
                         spec.v_div_vcrit)
    return isos


class _Job:
    def __init__(self, nml_file: str):
        self.nml_file = nml_file

    def __call__(self, input_file: str) -> bool:
        try:
            make_isochrones(input_file, self.nml_file)
        except MakeIsoError as exc:
            print(str(exc), file=sys.stderr)
            print(' problem in make_iso: no output written', file=sys.stderr)
            return False
        return True


def run(inputs: List[str], nml_file: str = 'input.nml',
        jobs: Optional[int] = 1, backend: str = 'processes') -> int:
    workers = resolve_workers(jobs, len(inputs))
    if workers > 1:
        print(f" make_iso: processing {len(inputs)} input files on "
              f"{workers} {backend}", file=sys.stderr)
    t0 = time.perf_counter()
    failed = 0
    for name, result in map_ordered(_Job(nml_file), inputs, workers, backend):
        if result is not True:
            failed += 1
            if isinstance(result, Exception):
                print(f"make_iso: problem processing {name}: {result}",
                      file=sys.stderr)
    if len(inputs) > 1:
        print(f" make_iso: {len(inputs) - failed}/{len(inputs)} grids in "
              f"{time.perf_counter() - t0:.1f} s", file=sys.stderr)
    return 1 if failed else 0


def main(argv: Optional[List[str]] = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if not argv:
        print('   make_iso                   ')
        print('   usage: ./make_iso.py <input> [<input> ...] [-j N]')
        print('   error: no command line argument provided')
        return 1
    parser = argparse.ArgumentParser(
        prog='make_iso.py',
        description='Build isochrones from EEP tracks (port of make_iso).')
    parser.add_argument('inputs', nargs='+',
                        help='input file(s); one isochrone set per file')
    parser.add_argument('--nml', default='input.nml',
                        help='namelist file holding &iso_controls '
                             '(default: ./input.nml, as in the Fortran)')
    add_parallel_args(parser, default=1)
    args = parser.parse_args(argv)
    return run(args.inputs, args.nml, args.jobs, args.backend)


if __name__ == '__main__':
    raise SystemExit(main())
