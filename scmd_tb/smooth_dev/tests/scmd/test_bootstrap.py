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


def test_fresh_reference():
    ctx = _ctx()
    a = B.realisations(ctx, 20001, 8.3, 13.28, 0.13, range(3), fresh_reference=True)
    b = B.realisations(ctx, 20001, 8.3, 13.28, 0.13, [2], fresh_reference=True)
    c = B.realisations(ctx, 20001, 8.3, 13.28, 0.13, range(3))
    assert a[2]["ks_statistic"] == b[0]["ks_statistic"]                 # depends on the number only
    assert [r["ks_statistic"] for r in a] != [r["ks_statistic"] for r in c]
    assert len({r["n_observed"] for r in a}) >= 1


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


def _ast_frame(n=6000, seed=0):
    rs = np.random.RandomState(seed)
    vin, iin = rs.uniform(-12, -9, n), rs.uniform(-12, -9, n)
    a = pd.DataFrame({"InputX": rs.uniform(2500, 3500, n), "InputY": rs.uniform(2500, 3500, n),
                      "InputF606W": vin, "InputF814W": iin})
    core = rs.normal(0, 0.02, n)
    blend = np.where(rs.random_sample(n) < 0.05, -rs.exponential(0.15, n), 0.0)   # a brightening tail
    a["OutputX"], a["OutputY"] = a["InputX"], a["InputY"]
    a["OutputF606W"] = vin + core + blend + rs.normal(0, 0.01, n)
    a["OutputF814W"] = iin + core + blend + rs.normal(0, 0.01, n)
    return a.astype(np.float32)


def test_smooth_draws():
    from isochronetoolbox.scmd.ast import NearestErrors
    a = _ast_frame()
    copy, smooth = NearestErrors(a, k=100), NearestErrors(a, k=100, smooth=True)
    q = np.random.RandomState(3)
    m = 20000
    args = (np.full(m, 300.0), q.uniform(-11.5, -9.5, m), q.uniform(-11.5, -9.5, m))
    ok_c, ev_c, ei_c = copy.draw(np.random.RandomState(4), *args)
    ok_s, ev_s, ei_s = smooth.draw(np.random.RandomState(4), *args)
    assert np.array_equal(ok_c, ok_s)                                      # same completeness draws
    assert len(np.unique(ev_s[ok_s])) == ok_s.sum()                         # continuous: no repeated values
    assert len(np.unique(ev_c[ok_c])) < 0.5 * ok_c.sum()                    # copies repeat
    mad = lambda x: 1.4826 * np.median(np.abs(x - np.median(x)))
    assert abs(mad(ev_s[ok_s]) / mad(ev_c[ok_c]) - 1) < 0.05                # core width kept
    assert abs(np.corrcoef(ev_s[ok_s], ei_s[ok_s])[0, 1] - np.corrcoef(ev_c[ok_c], ei_c[ok_c])[0, 1]) < 0.05
    assert abs(np.mean(ev_s[ok_s] > 0.1) - np.mean(ev_c[ok_c] > 0.1)) < 0.01   # blend tail kept (dV = in - out > 0)


def test_ast_resample():
    import os
    from isochronetoolbox.scmd.ast import ArtificialStars
    a = _ast_frame()
    fit = pd.DataFrame({"x": np.random.RandomState(1).uniform(2800, 3200, 500),
                        "y": np.random.RandomState(2).uniform(2800, 3200, 500)})
    zp = {"F606W": 30.0, "F814W": 29.0}
    win = {"F606W": (18.2, 20.8), "F814W": (17.2, 19.8)}
    full = ArtificialStars(a, fit, zp, win, {"F606W": 19.5, "F814W": 18.5}, k=50, smooth=True)
    r1, r2, r1b = full.resample(7), full.resample(8), full.resample(7)
    assert r1.n_ast == full.n_ast and np.array_equal(r1.nearest.ev, r1b.nearest.ev)
    assert not np.array_equal(r1.nearest.ev, r2.nearest.ev)
    assert len(np.unique(r1.nearest.ev)) < 0.7 * r1.n_ast                    # with replacement: ~37% repeats
    assert r1.nearest.smooth and r1.nearest.k == 50


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(name, "OK")
