"""Tests for isochronetoolbox.scmd.

Bit-for-bit agreement with the old three-script pipeline is checked by
``compare_with_legacy.py`` (it needs those scripts and the full isochrone
grid).  These tests need neither:

  - the cached KDE equals KDEpy's FFTKDE to the bit
  - one 3-band BC interpolator equals three single-band ones
  - the legacy 1 % age match and the exact one
  - a regenerated sCMD is the one that was scored
  - the no-survivor branch of the old stage 2, the worker's claiming and
    collect's completeness rule
  - a regression pin on a one-age fixture (exact numbers depend on the
    scipy build -- Qhull and pocketfft -- so it runs on scipy 1.13 only)

Run directly (``python tests/scmd/test_scmd.py``) or under pytest.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import scipy

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)

from isochronetoolbox.scmd import (  # noqa: E402
    KDE, BCTable, Context, Setup, process_age, regenerate, select_age)
from isochronetoolbox.scmd.pipeline import DATA  # noqa: E402

FIXTURE = os.path.join(HERE, "data", "iso_20001_8.3Gyr.txt")


def _setup() -> Setup:
    return Setup(iso=FIXTURE, dm_grid=(13.27, 13.315, 0.01), ebv_grid=(0.12, 0.165, 0.01))


def _kdepy_grid(cr, mr, n=1000, eps=1e-6):
    cm, mm = np.meshgrid(np.linspace(cr[0] - eps, cr[1] + eps, n), np.linspace(mr[0] - eps, mr[1] + eps, n))
    pts = np.column_stack([cm.ravel(), mm.ravel()])
    return pts[np.lexsort((pts[:, 1], pts[:, 0]))]


def test_kde_matches_kdepy():
    from KDEpy import FFTKDE
    rs = np.random.RandomState(3)
    for cr, mr in [((1.0, 1.5), (15.0, 19.0)), ((0.913, 1.688), (16.26, 17.99))]:
        kde, grid = KDE(cr, mr), _kdepy_grid(cr, mr)
        d = np.column_stack([rs.uniform(cr[0] + 0.01, cr[1] - 0.01, 400),
                             rs.uniform(mr[0] + 0.01, mr[1] - 0.01, 400)])
        ref = FFTKDE(kernel="gaussian", bw=0.15, norm=1).fit(d).evaluate(grid)
        assert np.array_equal(kde.pdf(d), ref)


def test_bc_one_interpolator_matches_three():
    from scipy.interpolate import LinearNDInterpolator
    path = os.path.join(DATA, "bc_mist_gaia_edr3_fehm300.npz")
    z = np.load(path)
    rs = np.random.RandomState(5)
    teff, logg = rs.uniform(3500, 7000, 300), rs.uniform(1.0, 5.0, 300)
    one = BCTable(path)._interp(teff, logg)
    for k in range(3):
        assert np.array_equal(one[:, k], LinearNDInterpolator(z["points"], z["bc"][:, k])(teff, logg))


def test_age_match():
    ages = np.repeat([9.8e9, 9.9e9, 10.0e9, 10.1e9, 10.2e9], 3)
    iso = pd.DataFrame({"isochrone_age_yr": ages, "initial_mass": np.arange(len(ages))})
    age, rows = select_age(iso, 10.0)
    assert age == 10.0e9 and sorted(set(rows["isochrone_age_yr"])) == [9.9e9, 10.0e9, 10.1e9]
    age, rows = select_age(iso, 10.0, exact=True)
    assert sorted(set(rows["isochrone_age_yr"])) == [10.0e9]
    age, rows = select_age(iso, 9.83)                  # closest file age, 1 % = 0.098 Gyr
    assert age == 9.8e9 and sorted(set(rows["isochrone_age_yr"])) == [9.8e9]


def test_regenerate_is_what_was_scored():
    ctx = Context(_setup())
    js, ks = process_age(ctx, 20001, 8.3)
    assert len(ks) == 4 and [r["kde_rank"] for r in ks] == [1, 2, 3, 4]
    assert len(js) == len(ctx.s.dm) * len(ctx.s.ebv)
    for r in ks[:2]:
        cat = regenerate(ctx, 20001, 8.3, r["distance_modulus"], r["extinction_ebv"])
        again = ctx.stage3(cat)
        assert again["ks_statistic"] == r["ks_statistic"]
        assert again["n_simulated"] == r["n_simulated"]
    js2, ks2 = process_age(ctx, 20001, 8.3)                     # deterministic
    assert [r["ks_statistic"] for r in ks2] == [r["ks_statistic"] for r in ks]


def test_no_survivor_branch():
    """No sCMD with 100 stars in the stage-2 window: every sCMD goes to stage 3."""
    s = _setup()
    s.kde_window = ((1.0, 1.002), (15.0, 19.0))           # far too narrow to hold 100 stars
    js, ks = process_age(Context(s), 20001, 8.3)
    assert js == [] and len(ks) == len(s.dm) * len(s.ebv)
    assert all(r["kde_rank"] == 0 and np.isnan(r["kde_js"]) for r in ks)


def test_worker_and_collect():
    """Unreadable isochrones are never claimed; collect takes complete runs only."""
    import tempfile
    from isochronetoolbox.scmd.pipeline import collect, status, worker
    with tempfile.TemporaryDirectory() as out:
        s = _setup()
        s.iso = os.path.join(HERE, "data", "iso_{run}_8.3Gyr.txt")
        s.ages, s.chunk = (8.3, 8.4), 1
        assert worker(s, out, [20001, 99999], log=lambda m: None) == 2
        st = status(s, out, [20001, 99999])
        assert len(st["done"]) == 2 and st["claimed"] == [] and len(st["todo"]) == 2
        os.remove(os.path.join(out, "units", "20001", "20001_8.4-8.4.ks.csv"))
        comb, rep = collect(s, out)
        assert rep["incomplete"] == [20001] and comb.empty


# filename, ks_statistic, js_divergence, n_simulated, kde_js of the four kept sCMDs
# (numpy 1.26.4, scipy 1.13.1, KDEpy 1.1.12)
PIN = [
    ("NGC6791_age8.3Gyr_DM13.270_EBV0.120.csv", 0.25818153587019316, 0.14243029832459542, 1115, 0.13023071722490287),
    ("NGC6791_age8.3Gyr_DM13.280_EBV0.120.csv", 0.278206256928138, 0.1414812417798668, 1117, 0.13561483678185388),
    ("NGC6791_age8.3Gyr_DM13.290_EBV0.120.csv", 0.26248834970259044, 0.15069135737115766, 1096, 0.13760840272677197),
    ("NGC6791_age8.3Gyr_DM13.270_EBV0.130.csv", 0.40598619480701814, 0.16297997425478863, 1119, 0.1428658553557794),
]


def test_regression_pin():
    if PIN is None or not scipy.__version__.startswith("1.13"):
        print("regression pin skipped (needs scipy 1.13)")
        return
    _, ks = process_age(Context(_setup()), 20001, 8.3)
    got = [(r["filename"], r["ks_statistic"], r["js_divergence"], r["n_simulated"], r["kde_js"]) for r in ks]
    for g, p in zip(got, PIN):
        assert g[0] == p[0] and g[3] == p[3], (g, p)
        assert np.allclose(g[1:3] + (g[4],), p[1:3] + (p[4],), rtol=1e-12, atol=0), (g, p)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(name, "OK")
