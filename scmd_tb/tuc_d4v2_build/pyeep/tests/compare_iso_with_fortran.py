#!/usr/bin/env python3
"""Regression: run the Fortran make_iso and this port on the same grids and
iso_controls settings, and compare the output files.

Needs a compiled Fortran ``make_iso`` (build it yourself from upstream ``iso``;
it is not distributed here) and existing ``.eep`` grids::

    python3 tests/compare_iso_with_fortran.py --fortran /path/to/make_iso \\
        --columns /path/to/my_history_columns_basic.list \\
        input_A.example input_B.example ...

Each input file is rerun under every configuration in CONFIGS, with output
redirected to a scratch directory.  A file is "identical" when the bytes
match; otherwise the largest relative difference over all numbers is shown.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
MAKE_ISO = os.path.join(os.path.dirname(HERE), 'make_iso.py')

# (name, &iso_controls body, ages block, optional trailing eep_style line)
LIN_AGES = 'min_max\nlinear\n81\n4.0e9\n12.0e9\n'
CONFIGS = [
    ('linear', 'do_smooth=.false.\ndo_PAV=.false.\ndo_linear_interpolation=.true.\n'
               'log_age_delta=0.2\nvery_low_mass_limit=0.6\n', LIN_AGES, ''),
    ('defaults', '', LIN_AGES, ''),
    ('steffen', 'do_smooth=.false.\ndo_PAV=.false.\n', LIN_AGES, ''),
    ('smooth_pav', 'do_smooth=.true.\ndo_PAV=.true.\nlog_age_delta=0.3\n',
     'list\nlog10\n5\n8.5\n9.0\n9.6\n9.9\n10.1\n', ''),
    ('max_eep', 'set_max_eep_number=.true.\nnew_max_eep_number=400\n'
                'do_linear_interpolation=.true.\n', LIN_AGES, ''),
    ('double', 'do_smooth=.false.\n', 'min_max\nlog10\n21\n9.0\n10.1\n', 'double'),
]


def rewrite_input(src: str, iso_dir: str, ages: str, style: str) -> str:
    with open(src) as fh:
        lines = fh.read().splitlines()
    lines[7] = iso_dir
    ntrk = int(lines[11].split()[0])
    head = lines[:12 + ntrk]
    tail = ['#specify isochrones', 'out.iso'] + ages.splitlines()
    if style:
        tail.append(style)
    return '\n'.join(head + tail) + '\n'


def numbers(path: str):
    out = []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith('#') or not ln.strip():
                continue
            out.extend(float(x) for x in ln.split())
    return out


def blocks(path: str):
    """Per isochrone: list of (eep, row text) data rows."""
    out, cur = [], None
    with open(path) as fh:
        for ln in fh:
            if ln.startswith('# number of EEPs'):
                cur = []
                out.append(cur)
            elif cur is not None and ln.strip() and not ln.startswith('#'):
                cur.append((int(ln[:5]), ln))
    return out


def dropped_only(a: str, b: str, dropped: set) -> bool:
    """True if ``b`` equals ``a`` except for rows at EEPs the port dropped.

    The header rows of each isochrone carry the EEP count, so only data rows
    are compared; every surviving row must match byte-for-byte.
    """
    ba, bb = blocks(a), blocks(b)
    if len(ba) != len(bb):
        return False
    for ra, rb in zip(ba, bb):
        if [r for r in ra if r[0] not in dropped] != [r for r in rb if r[0] not in dropped]:
            return False
        if len(rb) != len([r for r in ra if r[0] not in dropped]) + \
                len([r for r in rb if r[0] in dropped]):
            return False
    return True


def pav_pool_only(a: str, b: str, dropped: set) -> bool:
    """True if every difference sits in one run of rows touching a dropped EEP.

    With PAV on, the Fortran pools its stale value at a dropped EEP with its
    neighbours, so a contiguous run next to it differs as well.  (Feeding the
    Fortran's value into this port's PAV reproduces the Fortran rows exactly;
    see the module README.)
    """
    for ra, rb in zip(blocks(a), blocks(b)):
        rb_map = dict(rb)
        eeps = [e for e, _ in ra]
        bad = [i for i, (e, r) in enumerate(ra) if rb_map.get(e) != r]
        if not bad:
            continue
        if bad != list(range(bad[0], bad[-1] + 1)):
            return False
        if not ({eeps[bad[0]], eeps[bad[-1]]} & dropped):
            return False
        if set(rb_map) - set(eeps):
            return False
    return True


def compare(a: str, b: str, dropped=frozenset()) -> str:
    if not os.path.exists(a) and not os.path.exists(b):
        return 'both wrote nothing'
    if not (os.path.exists(a) and os.path.exists(b)):
        return f'MISMATCH: only one wrote output ({a if os.path.exists(a) else b})'
    with open(a, 'rb') as fa, open(b, 'rb') as fb:
        if fa.read() == fb.read():
            return 'identical'
    if dropped and dropped_only(a, b, dropped):
        return (f'identical except EEP(s) {sorted(dropped)} dropped for '
                f'repeated x values (Fortran reads stale memory there)')
    if dropped and pav_pool_only(a, b, dropped):
        return (f'identical except EEP(s) {sorted(dropped)} dropped and the '
                f'PAV pool next to it')
    na, nb = numbers(a), numbers(b)
    if len(na) != len(nb):
        return f'MISMATCH: {len(na)} vs {len(nb)} numbers'
    worst = max((abs(x - y) / max(abs(x), abs(y), 1e-300)
                 for x, y in zip(na, nb) if x != y), default=0.0)
    return f'DIFFERENT: worst relative difference {worst:.2e}'


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument('--fortran', required=True)
    ap.add_argument('--columns', required=True)
    ap.add_argument('--configs', default=','.join(c[0] for c in CONFIGS))
    ap.add_argument('inputs', nargs='+')
    args = ap.parse_args()
    wanted = args.configs.split(',')

    bad = 0
    known = 0
    for name, body, ages, style in (c for c in CONFIGS if c[0] in wanted):
        for src in args.inputs:
            with tempfile.TemporaryDirectory() as root:
                shutil.copy(args.columns, os.path.join(root, os.path.basename(args.columns)))
                with open(os.path.join(root, 'input.nml'), 'w') as fh:
                    fh.write('&iso_controls\n' + body + '/\n')
                res = {}
                for who in ('fortran', 'python'):
                    d = os.path.join(root, who)
                    os.makedirs(d)
                    inp = os.path.join(root, f'{who}.example')
                    with open(inp, 'w') as fh:
                        fh.write(rewrite_input(src, d, ages, style))
                    cmd = ([args.fortran, inp] if who == 'fortran'
                           else [sys.executable, MAKE_ISO, inp])
                    run = subprocess.run(cmd, cwd=root, capture_output=True, text=True)
                    res[who] = os.path.join(d, 'out.iso')
                    if who == 'python':
                        dropped = {int(m) for m in re.findall(
                            r'EEP (\d+) dropped: repeated x', run.stdout)}
                verdict = compare(res['fortran'], res['python'], dropped)
                bad += not verdict.startswith(('identical', 'both'))
                known += verdict.startswith('identical except')
                print(f'{name:11s} {os.path.basename(src):40s} {verdict}')
    print(('ALL MATCH' if bad == 0 else f'{bad} MISMATCHES')
          + (f' ({known} with the documented repeated-x divergence)' if known else ''))
    return 1 if bad else 0


if __name__ == '__main__':
    raise SystemExit(main())
