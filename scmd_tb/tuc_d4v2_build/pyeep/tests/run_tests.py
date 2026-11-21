#!/usr/bin/env python3
"""Run the whole pyeep self-check suite and print a pass/fail summary.

    python3 tests/run_tests.py            # everything (a few seconds)
    python3 tests/run_tests.py iso        # only tests whose name contains 'iso'

Needs only NumPy -- no pytest, MESA data or Fortran.  Every ``test_*``
function in every ``tests/test_*.py`` file is run; the exit status is non-zero
if any fails, so this also serves as the CI entry point.
"""

from __future__ import annotations

import glob
import importlib.util
import inspect
import io
import os
import sys
import time
import traceback
from contextlib import redirect_stderr, redirect_stdout

HERE = os.path.dirname(os.path.abspath(__file__))


def main(argv) -> int:
    pattern = argv[0] if argv else ''
    sys.path.insert(0, HERE)
    sys.path.insert(0, os.path.dirname(HERE))
    passed, failed = 0, []
    t_all = time.perf_counter()
    for path in sorted(glob.glob(os.path.join(HERE, 'test_*.py'))):
        name = os.path.splitext(os.path.basename(path))[0]
        spec = importlib.util.spec_from_file_location(name, path)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod   # so worker processes can unpickle its functions
        spec.loader.exec_module(mod)
        # only test functions defined in this file (not imported helpers)
        tests = sorted(n for n, obj in vars(mod).items()
                       if n.startswith('test_') and inspect.isfunction(obj)
                       and obj.__module__ == name)
        for fn in tests:
            label = f"{name}::{fn}"
            if pattern not in label:
                continue
            log = io.StringIO()
            t0 = time.perf_counter()
            try:
                with redirect_stdout(log), redirect_stderr(log):
                    getattr(mod, fn)()
            except Exception:  # noqa: BLE001 - report and keep going
                failed.append((label, log.getvalue() + traceback.format_exc()))
                status = 'FAIL'
            else:
                passed += 1
                status = 'pass'
            print(f"  {status}  {label}  ({time.perf_counter() - t0:.2f} s)")

    for label, detail in failed:
        print(f"\n---- {label} ----\n{detail}")
    total = passed + len(failed)
    print(f"\n{passed}/{total} passed in {time.perf_counter() - t_all:.1f} s")
    return 1 if failed or total == 0 else 0


if __name__ == '__main__':
    raise SystemExit(main(sys.argv[1:]))
