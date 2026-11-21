#!/usr/bin/env python3
"""Regression-test the Python port against the Fortran ``make_eep`` binary.

Runs both programs on the *same* two-part input file, writing to two scratch
directories, and compares the resulting ``.eep`` files header-by-header and
number-by-number.

Usage::

    cd <iso repo root>          # where input.nml / input.eep / input.format live
    python3 python/tests/compare_with_fortran.py input.example

Options::

    --fortran PATH   path to the Fortran make_eep binary (default ./make_eep)
    --workdir PATH   scratch directory (default: a temporary directory)
    --rtol FLOAT     relative tolerance for the numeric comparison (default 1e-15)
    --keep           do not delete the scratch directory on exit

Exit status is 0 when every track matches, 1 otherwise.

Note on ``initial_mass``
------------------------
The Fortran ``track`` derived type leaves ``initial_mass`` uninitialised, and
``set_star_type_from_history`` is called before the value is recovered from the
track data.  When a MESA history header carries no ``initial_mass`` field the
Fortran therefore tests recycled heap memory (in practice the *previous*
track's mass, or garbage for the first track).  This only changes the outcome
for tracks whose central hydrogen varies by less than 0.1 over the whole track,
where the sub-stellar branch of ``ZAMS`` may be taken spuriously.  The Python
port resolves ``initial_mass`` from the data first, so it is deterministic and
order-independent; expect a mismatch against the Fortran for such tracks unless
they are preceded in the input list by a track of mass >= 0.1.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

REPO_PY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def rewrite_input(src: str, dst: str, eep_dir: str) -> None:
    """Copy the two-part input file, redirecting the eeps directory."""
    with open(src) as fh:
        lines = fh.read().split('\n')
    # line 4 (0-based) is the '#data directories' comment; 5,6,7 are the
    # history / eeps / isochrones paths
    lines[6] = eep_dir
    with open(dst, 'w') as fh:
        fh.write('\n'.join(lines))


def read_eep(path: str):
    header, rows = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith('#'):
                header.append(line.rstrip())
            elif line.strip():
                rows.append([float(x) for x in line.split()])
    return header, np.array(rows)


def compare(fdir: str, pdir: str, rtol: float) -> int:
    fnames = sorted(f for f in os.listdir(fdir) if f.endswith('.eep'))
    pnames = sorted(f for f in os.listdir(pdir) if f.endswith('.eep'))
    if fnames != pnames:
        only_f = sorted(set(fnames) - set(pnames))
        only_p = sorted(set(pnames) - set(fnames))
        print(f"FAIL: different output files (fortran-only={only_f} python-only={only_p})")
        return 1

    failures = 0
    worst = 0.0
    for name in fnames:
        hf, af = read_eep(os.path.join(fdir, name))
        hp, ap = read_eep(os.path.join(pdir, name))

        # the Fortran pads the column-name rows with trailing blanks
        hdr_ok = [h.rstrip() for h in hf] == [h.rstrip() for h in hp]
        if not hdr_ok:
            print(f"FAIL {name}: header mismatch")
            for a, b in zip(hf, hp):
                if a.rstrip() != b.rstrip():
                    print(f"      fortran: {a.rstrip()}")
                    print(f"      python : {b.rstrip()}")
            failures += 1
            continue

        if af.shape != ap.shape:
            print(f"FAIL {name}: shape {af.shape} (fortran) vs {ap.shape} (python)")
            failures += 1
            continue

        absdiff = np.abs(af - ap)
        rel = np.where(np.abs(af) > 0, absdiff / np.maximum(np.abs(af), 1e-300), absdiff)
        maxrel = float(rel.max()) if rel.size else 0.0
        worst = max(worst, maxrel)
        if maxrel > rtol:
            i, j = np.unravel_index(np.argmax(rel), rel.shape)
            print(f"FAIL {name}: max rel diff {maxrel:.3e} > {rtol:.1e} "
                  f"at row {i}, col {j} ({af[i, j]!r} vs {ap[i, j]!r})")
            failures += 1
        else:
            print(f"OK   {name}: rows={af.shape[0]} max rel diff {maxrel:.3e}")

    print()
    print(f"{len(fnames) - failures}/{len(fnames)} tracks match "
          f"(worst relative difference {worst:.3e}, tolerance {rtol:.1e})")
    return 1 if failures else 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input', help='two-part input file (e.g. input.example)')
    ap.add_argument('--fortran', default='./make_eep', help='Fortran make_eep binary')
    ap.add_argument('--workdir', default=None, help='scratch directory')
    ap.add_argument('--rtol', type=float, default=1e-15, help='relative tolerance')
    ap.add_argument('--keep', action='store_true', help='keep the scratch directory')
    ap.add_argument('-j', '--jobs', type=int, default=1, metavar='N',
                    help='parallel workers for the Python run; 0 = auto-detect '
                         '(the Fortran binary is always serial)')
    args = ap.parse_args(argv)

    if not os.path.exists(args.fortran):
        print(f"error: Fortran binary not found: {args.fortran}\n"
              f"       build it with 'cd make && make' (needs MESA_DIR), or pass --fortran",
              file=sys.stderr)
        return 2

    work = args.workdir or tempfile.mkdtemp(prefix='eepcmp_')
    fdir = os.path.join(work, 'fortran')
    pdir = os.path.join(work, 'python')
    for d in (fdir, pdir):
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)

    fin = os.path.join(work, 'in_fortran.example')
    pin = os.path.join(work, 'in_python.example')
    rewrite_input(args.input, fin, fdir)
    rewrite_input(args.input, pin, pdir)

    try:
        print(f"running {args.fortran} ...")
        rc = subprocess.run([os.path.abspath(args.fortran), fin],
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode
        if rc != 0:
            print(f"error: Fortran make_eep exited {rc}", file=sys.stderr)
            return 2

        print("running python make_eep ...")
        cmd = [sys.executable, os.path.join(REPO_PY, 'make_eep.py'), pin]
        if args.jobs != 1:
            cmd += ['-j', str(args.jobs)]
        rc = subprocess.run(cmd,
                            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode
        if rc != 0:
            print(f"error: python make_eep exited {rc}", file=sys.stderr)
            return 2

        print()
        return compare(fdir, pdir, args.rtol)
    finally:
        if args.keep:
            print(f"\nscratch directory kept at {work}")
        elif args.workdir is None:
            shutil.rmtree(work, ignore_errors=True)


if __name__ == '__main__':
    raise SystemExit(main())
