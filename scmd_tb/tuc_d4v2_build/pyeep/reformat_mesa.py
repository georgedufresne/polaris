#!/usr/bin/env python3
"""Rewrite MCDSEP-style track files into the MESA history format make_eep reads.

MCDSEP files store one value per line: after a six-line header the file is a
flat stream of numbers, one column-value per line, repeating for every model.
``make_eep`` expects standard MESA layout, one model per row.  This script
converts the former into the latter.

Usage::

    ./reformat_mesa.py <file-or-directory> [...] [-o OUTDIR]

Each input may be a file or a directory (all regular files inside are
converted, non-recursively).  By default the output is written alongside the
input with a ``_mesa23`` suffix; pass ``-o`` to collect results elsewhere.

Examples::

    # convert one file in place (-> m0700.feh031.20001_mesa23)
    ./reformat_mesa.py m0700.feh031.20001

    # convert a whole extracted directory into a clean output tree
    ./reformat_mesa.py /scratch/George/Files_For_Iso_Construction/TEMP20001 \
        -o /scratch/George/tracks_mesa

    # then point make_eep at it
    ./make_eep.py input.example
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pyeep.reformat import (  # noqa: E402
    DEFAULT_COLUMNS,
    is_mesa_format,
    reformat_mesa_file,
)
from pyeep.parallel import (  # noqa: E402
    add_parallel_args,
    map_ordered,
    resolve_workers,
)

DEFAULT_SUFFIX = '_mesa23'


def gather(paths, suffix):
    """Expand the command-line paths into a sorted list of input files."""
    files = []
    for p in paths:
        if os.path.isdir(p):
            for name in sorted(os.listdir(p)):
                full = os.path.join(p, name)
                if os.path.isfile(full) and not name.endswith(suffix):
                    files.append(full)
        elif os.path.isfile(p):
            files.append(p)
        else:
            print(f"warning: no such file or directory: {p}", file=sys.stderr)
    return files


class _ReformatJob:
    """One file's conversion, in a picklable form for a worker pool."""

    def __init__(self, outdir, suffix, columns, encoding, force, quiet):
        self.outdir = outdir
        self.suffix = suffix
        self.columns = columns
        self.encoding = encoding
        self.force = force
        self.quiet = quiet

    def __call__(self, src):
        if not self.force and is_mesa_format(src, self.encoding):
            if not self.quiet:
                print(f"skip {os.path.basename(src)}: already MESA format")
            return 'skipped'

        if self.outdir:
            dst = os.path.join(self.outdir, os.path.basename(src))
        else:
            dst = src + self.suffix

        nrow = reformat_mesa_file(src, dst, columns=self.columns,
                                  encoding=self.encoding)
        if not self.quiet:
            print(f"{os.path.basename(src)} -> {dst}  ({nrow} models)")
        return 'converted'


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paths', nargs='+', help='input files and/or directories')
    ap.add_argument('-o', '--outdir', default=None,
                    help='write results here instead of beside the input')
    ap.add_argument('--suffix', default=DEFAULT_SUFFIX,
                    help=f'suffix for in-place output (default {DEFAULT_SUFFIX})')
    ap.add_argument('--columns', default=None,
                    help='comma-separated column names (default: read from the '
                         'file header, else the standard 22-column layout)')
    ap.add_argument('--encoding', default=None,
                    help='force a text encoding instead of detecting one')
    ap.add_argument('--force', action='store_true',
                    help='convert even files that already look like MESA format')
    ap.add_argument('-q', '--quiet', action='store_true',
                    help='only report errors and the final summary')
    add_parallel_args(ap, default=1)
    args = ap.parse_args(argv)

    columns = None
    if args.columns:
        columns = [c.strip() for c in args.columns.split(',') if c.strip()]
        if not columns:
            print("error: --columns is empty", file=sys.stderr)
            return 2

    files = gather(args.paths, args.suffix)
    if not files:
        print("error: no input files found", file=sys.stderr)
        return 2

    if args.outdir:
        os.makedirs(args.outdir, exist_ok=True)

    workers = resolve_workers(args.jobs, len(files))
    if workers > 1 and not args.quiet:
        print(f"converting {len(files)} files on {workers} {args.backend}",
              file=sys.stderr)

    job = _ReformatJob(args.outdir, args.suffix, columns, args.encoding,
                       args.force, args.quiet)

    done = skipped = failed = 0
    for src, result in map_ordered(job, files, workers, args.backend):
        if isinstance(result, Exception):
            print(f"FAIL {os.path.basename(src)}: {result}", file=sys.stderr)
            failed += 1
        elif result == 'skipped':
            skipped += 1
        else:
            done += 1

    print(f"\nconverted {done}, skipped {skipped}, failed {failed}")
    return 1 if failed else 0


if __name__ == '__main__':
    raise SystemExit(main())
