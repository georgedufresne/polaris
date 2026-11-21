"""Tests for the artificial-star error model options (ast.py) and common random numbers.

  - metric="crowding" keeps short- and long-exposure artificial stars apart
  - radial="corrected" undoes radius-dependent incompleteness of the fit stars
  - common random numbers: every grid cell starts from the same seed, regenerate
    returns the scored catalogue

Synthetic inputs only.  Run directly (``python tests/scmd/test_ast_model.py``) or under pytest.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from isochronetoolbox.scmd import Context, Setup, process_age, regenerate  # noqa: E402
from isochronetoolbox.scmd.ast import ArtificialStars, NearestErrors       # noqa: E402

ZP = {"F606W": 30.0, "F814W": 29.0}


def _asts(n=20000, seed=0, comp=None):
    """Artificial stars over V_vega 15.5-19.5, radius 0-2500 px around (3000, 3000);
    ``comp(r)`` sets the recovery probability."""
    rs = np.random.RandomState(seed)
    r = 2500 * np.sqrt(rs.random_sample(n))
    t = rs.uniform(0, 2 * np.pi, n)
    x, y = 3000 + r * np.cos(t), 3000 + r * np.sin(t)
    v = rs.uniform(15.5, 19.5, n) - ZP["F606W"]
    i = v + ZP["F606W"] - 0.6 - ZP["F814W"]
    e = rs.normal(0, 0.01 + 0.03 * (r < 500), n)
    lost = rs.random_sample(n) > (comp(r) if comp else 1.0)
    out = pd.DataFrame({"InputX": x, "InputY": y, "InputF606W": v, "InputF814W": i,
                        "OutputX": x + 2.0 * lost, "OutputY": y,
                        "OutputF606W": v - e, "OutputF814W": i - e})
    return out.astype(np.float32)


def test_crowding_metric_keeps_exposures_apart():
    a = _asts()
    split = (16.48 - ZP["F606W"], 15.55 - ZP["F814W"])
    ne = NearestErrors(a, k=30, metric="crowding", split=split)
    vin = a["InputF606W"].values
    short = vin < split[0]
    mixed = short[ne.nbr].any(axis=1) & (~short[ne.nbr]).any(axis=1)
    assert not mixed.any()
    lin = NearestErrors(a, k=30)
    assert (short[lin.nbr].any(axis=1) & (~short[lin.nbr]).any(axis=1)).any()   # the linear metric mixes them


def test_corrected_radii_undo_incompleteness():
    comp = lambda r: 0.5 + 0.5 * np.clip(r / 1500.0, 0, 1)        # half the core stars lost
    a = _asts(60000, comp=comp)
    rs = np.random.RandomState(1)
    n = 40000
    r_true = 2500 * np.sqrt(rs.random_sample(n))                    # uniform surface density
    t = rs.uniform(0, 2 * np.pi, n)
    keep = rs.random_sample(n) < comp(r_true)                       # the observed (incomplete) fit stars
    v = rs.uniform(16.6, 19.4, n)[keep]
    fit = pd.DataFrame({"x": 3000 + (r_true * np.cos(t))[keep], "y": 3000 + (r_true * np.sin(t))[keep],
                        "v": v, "i": v - 0.6})
    win = {"F606W": (15.5, 19.5), "F814W": (14.9, 18.9)}
    sl = {"F606W": 16.48, "F814W": 15.55}
    obs = ArtificialStars(a, fit, ZP, win, sl, k=30, metric="crowding")
    cor = ArtificialStars(a, fit, ZP, win, sl, k=30, metric="crowding", radial="corrected")
    draw = lambda m: m.radii(np.random.RandomState(2), 50000, np.full(50000, 18.0))
    frac_core = lambda r: np.mean(r < 750)
    truth = (750 / 2500) ** 2
    assert abs(frac_core(draw(cor)) - truth) < 0.01                 # corrected: the true core fraction
    assert frac_core(draw(obs)) < truth - 0.02                       # observed: core depleted


def _ctx(crn):
    s = Setup(iso=os.path.join(HERE, "data", "iso_20001_8.3Gyr.txt"), dm_grid=(13.27, 13.315, 0.01),
              ebv_grid=(0.12, 0.145, 0.01), crn=crn)
    return Context(s)


def test_common_random_numbers():
    ctx = _ctx(True)
    _, make = ctx.model(20001, 8.3)
    cells = list(ctx.grid(20001, 8.3)[1])
    for dm, ebv, cat in cells[:3]:
        assert np.array_equal(cat, make(np.random.RandomState(ctx.s.crn_seed), dm, ebv))
        assert np.array_equal(cat, regenerate(ctx, 20001, 8.3, dm, ebv))
    js, ks = process_age(ctx, 20001, 8.3)
    r = ks[0]
    again = ctx.stage3(regenerate(ctx, 20001, 8.3, r["distance_modulus"], r["extinction_ebv"]))
    assert again["ks_statistic"] == r["ks_statistic"]
    old = _ctx(False)                                                 # without crn: the legacy stream
    assert not np.array_equal(list(old.grid(20001, 8.3)[1])[1][2], cells[1][2])


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(name, "OK")
