"""Shared helpers for running pyeep work across multiple CPUs.

Why processes and not threads
-----------------------------
The dominant cost in ``make_eep`` is parsing history files — ``str.split`` and
``float()`` over millions of fields. That is pure-Python work holding the
GIL, so threads do not overlap; measured on a 16-core Ryzen 9 7950X with eight
large tracks:

===========  ========
backend      wall
===========  ========
serial       5.2 s
threads      7.4 s  (slower than serial: GIL contention plus overhead)
processes    0.8 s  (6.5x speedup)
===========  ========

Processes are therefore the default. A ``threads`` backend is still selectable
via ``--backend`` so the above can be reproduced, but it is not recommended.

Each track is fully independent — its own input file, its own output file, no
shared mutable state — so the work is embarrassingly parallel and results are
identical regardless of worker count.
"""

from __future__ import annotations

import argparse
import io
import os
import sys
from contextlib import redirect_stdout, redirect_stderr
from typing import Callable, Iterable, Iterator, List, Optional, Sequence, Tuple

#: Conservative peak resident memory per worker, in MiB.  Measured peak for a
#: ~90 MB / 100k-model history file is roughly 190 MiB; round up for headroom.
MEM_PER_WORKER_MB = 300


def available_memory_mb() -> Optional[int]:
    """Return MemAvailable in MiB, or ``None`` if it cannot be determined."""
    try:
        with open('/proc/meminfo') as fh:
            for line in fh:
                if line.startswith('MemAvailable:'):
                    return int(line.split()[1]) // 1024
    except (OSError, ValueError, IndexError):
        pass
    return None


def cpu_count() -> int:
    """Usable CPU count, honouring CPU affinity / cgroup limits when visible."""
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except AttributeError:  # pragma: no cover - non-Linux
        return max(1, os.cpu_count() or 1)


def resolve_workers(requested: Optional[int],
                    n_items: int,
                    mem_per_worker_mb: int = MEM_PER_WORKER_MB,
                    quiet: bool = False) -> int:
    """Decide how many workers to use.

    ``requested`` is the raw ``--jobs`` value: ``None`` or ``0`` means "choose
    automatically", a negative value means "all CPUs minus that many".  The
    result is clamped to the number of work items, to the CPU count, and — for
    the automatic case — to what available memory can support.
    """
    ncpu = cpu_count()

    if requested is None or requested == 0:
        workers = ncpu
        mem = available_memory_mb()
        if mem is not None:
            # leave ~1 GiB for the parent process and the page cache
            mem_workers = max(1, (mem - 1024) // mem_per_worker_mb)
            if mem_workers < workers:
                if not quiet:
                    print(f"note: limiting to {mem_workers} workers "
                          f"({mem} MiB available, ~{mem_per_worker_mb} MiB per worker); "
                          f"override with -j N",
                          file=sys.stderr)
                workers = mem_workers
    elif requested < 0:
        workers = max(1, ncpu + requested)
    else:
        workers = requested

    return max(1, min(workers, n_items, ncpu))


def add_parallel_args(parser: argparse.ArgumentParser,
                      default: Optional[int] = 1) -> None:
    """Add ``-j/--jobs`` and ``--backend`` to an argparse parser.

    ``default`` is the value used when the flag is absent; 1 keeps the
    single-process behaviour that matches the Fortran driver exactly.
    """
    parser.add_argument(
        '-j', '--jobs', type=int, default=default, metavar='N',
        help="number of parallel worker processes; 0 = auto-detect "
             "(CPU count, capped by available memory), negative = all CPUs "
             f"minus N (default: {default})")
    parser.add_argument(
        '--backend', choices=('processes', 'threads'), default='processes',
        help="parallel backend (default: processes; threads are slower here "
             "because history parsing is GIL-bound, and are provided only for "
             "comparison)")


class _Captured:
    """A worker result plus whatever the worker printed."""

    __slots__ = ('value', 'out', 'err', 'error')

    def __init__(self, value=None, out='', err='', error=None):
        self.value = value
        self.out = out
        self.err = err
        self.error = error


def _run_captured(payload):
    """Execute ``func(item)`` in a worker, capturing its stdout/stderr.

    Output is buffered rather than written directly so the parent can replay it
    in input order; otherwise concurrent workers would interleave their log
    lines and the run would be non-reproducible.
    """
    func, item = payload
    out, err = io.StringIO(), io.StringIO()
    try:
        with redirect_stdout(out), redirect_stderr(err):
            value = func(item)
        return _Captured(value, out.getvalue(), err.getvalue())
    except Exception as exc:  # noqa: BLE001 - reported by the parent
        return _Captured(None, out.getvalue(), err.getvalue(), exc)


def map_ordered(func: Callable,
                items: Sequence,
                workers: int,
                backend: str = 'processes') -> Iterator[Tuple[object, object]]:
    """Apply ``func`` to ``items``, yielding ``(item, result)`` in input order.

    Worker output is replayed on the parent's streams as each result is
    consumed, so the log is identical to a serial run regardless of the number
    of workers.  An exception inside a worker is yielded as the result rather
    than raised, letting the caller decide whether to continue the batch.

    With ``workers == 1`` the work runs inline — no pool, no pickling — so the
    default path stays byte-for-byte equivalent to the original serial driver.
    """
    if workers <= 1:
        for item in items:
            try:
                yield item, func(item)
            except Exception as exc:  # noqa: BLE001
                yield item, exc
        return

    if backend == 'threads':
        from concurrent.futures import ThreadPoolExecutor as Pool
    else:
        from concurrent.futures import ProcessPoolExecutor as Pool

    payloads = [(func, item) for item in items]
    with Pool(max_workers=workers) as pool:
        # executor.map preserves input order
        for item, cap in zip(items, pool.map(_run_captured, payloads)):
            if cap.out:
                sys.stdout.write(cap.out)
            if cap.err:
                sys.stderr.write(cap.err)
            yield item, (cap.error if cap.error is not None else cap.value)
