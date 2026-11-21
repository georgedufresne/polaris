"""Tests for the bootstrap (bootstrap.py) and the nearest-artificial-star errors (ast.py).

  - a realisation depends only on its bootstrap number: any split into slices
    gives the same rows, and the summary is mean + 3 x std (ddof = 0)
  - the slices of one directory must agree on the cell
  - NearestErrors returns one artificial star's own pair of errors

Run directly (``python tests/scmd/test_bootstrap.py``) or under pytest.
"""

from __future__ import annotations

import os
import sys
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from isochronetoolbox.scmd import Context, Setup  # noqa: E402
from isochronetoolbox.scmd import bootstrap as B  # noqa: E402

FIXTURE = os.path.join(HERE, "data", "iso_20001_8.3Gyr.txt")
CELL = {"simulation_id": 20001, "age_gyr": 8.3, "distance_modulus": 13.28, "extinction_ebv": 0.13}


def _ctx():
    return Context(Setup(iso=FIXTURE, dm_grid=(13.27, 13.315, 0.01), ebv_grid=(0.12, 0.165, 0.01)))


def test_slices_are_independent():
    ctx = _ctx()
    one = B.realisations(ctx, 20001, 8.3, 13.28, 0.13, range(6))
    with tempfile.TemporaryDirectory() as out:
        for task in range(3):
            B.run_slice(ctx, out, CELL, 6, task, 3)
        rep = B.collect(out)
        got = pd.read_csv(os.path.join(out, "bootstrap_2d_ks_results.csv"))
    assert rep["missing"] == 0 and list(got["bootstrap_number"]) == list(range(6))
    assert np.allclose(got["ks_statistic"].values, [r["ks_statistic"] for r in one], rtol=0, atol=1e-15)  # CSV round trip
    ks = got["ks_statistic"].values
    assert abs(rep["ks_3sigma"] - (ks.mean() + 3 * ks.std())) < 1e-15
    assert len({r["ks_statistic"] for r in one}) == 6            # the realisations differ


def test_cell_lock():
    ctx = _ctx()
    with tempfile.TemporaryDirectory() as out:
        B.run_slice(ctx, out, CELL, 2, 0, 2)
        try:
            B.run_slice(ctx, out, {**CELL, "distance_modulus": 13.29}, 2, 1, 2)
        except SystemExit:
            return
    raise AssertionError("a second cell was accepted")


def test_nearest_errors_pairs():
    from isochronetoolbox.scmd.ast import NearestErrors
    rs = np.random.RandomState(0)
    n = 4000
    vin, iin = rs.uniform(-12, -9, n), rs.uniform(-12, -9, n)
    a = pd.DataFrame({"InputX": rs.uniform(2500, 3500, n), "InputY": rs.uniform(2500, 3500, n),
                      "InputF606W": vin, "InputF814W": iin})
    shift = rs.normal(0, 0.02, n)
    a["OutputX"], a["OutputY"] = a["InputX"], a["InputY"]
    a["OutputF606W"], a["OutputF814W"] = vin + shift, iin + shift            # perfectly correlated errors
    a.loc[:99, "OutputF606W"] += 2.0                                       # 100 not recovered
    ne = NearestErrors(a.astype(np.float32), k=10)
    ok, ev, ei = ne.draw(rs, np.full(3000, 300.0), rs.uniform(-11.5, -9.5, 3000), rs.uniform(-11.5, -9.5, 3000))
    assert np.allclose(ev[ok], ei[ok], atol=1e-6)                          # both bands from one star
    assert 0.9 < ok.mean() < 1.0


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(name, "OK")
