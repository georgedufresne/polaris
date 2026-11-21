"""Tests for the parallel helpers (``pyeep.parallel``).

Checks worker-count resolution and that ``map_ordered`` preserves input order,
replays worker output deterministically, and surfaces worker exceptions rather
than raising or aborting the batch.

Run directly (``python tests/test_parallel.py``) or under pytest.
"""

from __future__ import annotations

import io
import os
import sys
from contextlib import redirect_stdout

HERE = os.path.dirname(os.path.abspath(__file__))
PY_ROOT = os.path.dirname(HERE)
sys.path.insert(0, PY_ROOT)
sys.path.insert(0, HERE)

from pyeep.parallel import (  # noqa: E402
    cpu_count,
    map_ordered,
    resolve_workers,
)


# module-level so the process backend can pickle them

def _square(x):
    return x * x


def _noisy(x):
    print(f"item {x}")
    return x


def _boom(x):
    if x == 3:
        raise ValueError("boom at 3")
    return x


def _pid_and_square(x):
    return os.getpid(), x * x


def test_resolve_workers():
    ncpu = cpu_count()

    # explicit request, clamped by item count and CPU count
    assert resolve_workers(4, 100) == min(4, ncpu)
    assert resolve_workers(100, 3) == min(3, ncpu)  # runners may have < 3 CPUs
    assert resolve_workers(1, 100) == 1

    # never returns less than 1
    assert resolve_workers(-10 ** 6, 100) == 1
    assert resolve_workers(5, 0) == 1

    # negative = all CPUs minus N
    assert resolve_workers(-1, 100) == max(1, ncpu - 1)

    # auto never exceeds CPU count or item count
    auto = resolve_workers(0, 100, quiet=True)
    assert 1 <= auto <= ncpu, auto
    assert resolve_workers(None, 2, quiet=True) <= 2

    # a tiny memory budget forces a single worker
    assert resolve_workers(0, 100, mem_per_worker_mb=10 ** 9, quiet=True) == 1
    print("resolve_workers OK")


def _check_order(backend, workers):
    items = list(range(25))
    got = list(map_ordered(_square, items, workers, backend))
    assert [i for i, _ in got] == items, "input order not preserved"
    assert [r for _, r in got] == [i * i for i in items], "wrong results"


def test_order_preserved():
    for backend in ('processes', 'threads'):
        for workers in (1, 4):
            _check_order(backend, workers)
    print("order preserved OK")


def test_output_replayed_in_order():
    """Worker prints must appear in input order, not completion order."""
    items = list(range(20))
    expected = ''.join(f"item {i}\n" for i in items)
    for backend in ('processes', 'threads'):
        buf = io.StringIO()
        with redirect_stdout(buf):
            list(map_ordered(_noisy, items, 4, backend))
        assert buf.getvalue() == expected, (
            f"{backend}: output out of order\n{buf.getvalue()!r}")
    print("output ordering OK")


def test_worker_exception_surfaced():
    """A failing item yields its exception; the rest of the batch still runs."""
    items = list(range(6))
    for backend in ('processes', 'threads'):
        for workers in (1, 3):
            got = dict(map_ordered(_boom, items, workers, backend))
            assert isinstance(got[3], ValueError), (backend, workers, got[3])
            for i in items:
                if i != 3:
                    assert got[i] == i, (backend, workers, i, got[i])
    print("worker exceptions OK")


def test_serial_path_runs_inline():
    """workers == 1 must not spawn a pool (keeps the default path unchanged)."""
    items = [1, 2, 3]
    pid = os.getpid()
    got = list(map_ordered(_pid_and_square, items, 1))
    assert [r[1] for _, r in got] == [1, 4, 9]
    # every item ran in this very process, so no pool was created
    assert all(r[0] == pid for _, r in got)
    print("serial inline OK")


if __name__ == '__main__':
    test_resolve_workers()
    test_order_preserved()
    test_output_replayed_in_order()
    test_worker_exception_surfaced()
    test_serial_path_runs_inline()
    print("\nall parallel tests passed")
