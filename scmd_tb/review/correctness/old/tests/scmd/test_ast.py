"""Tests for the artificial-star sCMD model (isochronetoolbox.scmd.ast).

  - lininterp equals a line-by-line transcription of MakeTestCMD's
    lininterp/locate, including on a non-monotonic array, except where the
    Fortran reads outside the array (between the first two points of a
    descending array) -- there it takes the first segment
  - read_dsep drops a block cut short, as in the truncated 47Tuc_iso files
  - the error draw honours completeness and stays inside each bin's errors
  - the binning reproduces M. Ying's 47 Tuc error files bin for bin (runs only
    where his data are present, /dartfs-hpc/scratch/f007058/Files_from_Martin)

Run directly (``python tests/scmd/test_ast.py``) or under pytest.
"""

from __future__ import annotations

import os
import sys
import tempfile

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from isochronetoolbox.scmd.ast import ArtificialStars, BandErrors, lininterp  # noqa: E402
from isochronetoolbox.scmd.model import read_dsep  # noqa: E402

MARTIN = "/dartfs-hpc/scratch/f007058/Files_from_Martin"


def _fortran_lininterp(x, y, xv):
    n, asc = len(x), x[-1] >= x[0]

    def locate(t):
        jl, ju = 0, n + 1
        while ju - jl > 1:
            jm = (ju + jl) // 2
            if asc == (t >= x[jm - 1]):
                jl = jm
            else:
                ju = jm
        return 1 if t == x[0] else (n - 1 if t == x[-1] else jl)

    out, oob = [], []
    for t in xv:
        if asc:
            if t <= x[0]:
                out.append(y[0]); oob.append(False); continue
            if t >= x[-1]:
                out.append(y[-1]); oob.append(False); continue
            i = locate(t)
            out.append((y[i] - y[i - 1]) / (x[i] - x[i - 1]) * (t - x[i - 1]) + y[i - 1]); oob.append(False)
        else:
            if t <= x[-1]:
                out.append(y[-1]); oob.append(False); continue
            if t >= x[0]:
                out.append(y[0]); oob.append(False); continue
            i = locate(t)
            oob.append(i == 1)
            out.append((y[i - 1] - y[i - 2]) / (x[i - 1] - x[i - 2]) * (t - x[i - 2]) + y[i - 2])
    return np.array(out), np.array(oob)


def test_lininterp_matches_fortran():
    rs = np.random.RandomState(1)
    x = np.array([18.5, 18.0, 17.6, 17.3, 17.36, 17.1, 16.5, 15.8, 15.0])   # descending, a wiggle
    y = np.linspace(0.6, 0.75, len(x))
    xv = rs.uniform(14.5, 19.0, 3000)
    ref, oob = _fortran_lininterp(x, y, xv)
    assert np.array_equal(lininterp(x, y, xv)[~oob], ref[~oob]) and oob.any()
    xa = np.sort(rs.uniform(0, 5, 25))
    xva = rs.uniform(-1, 6, 3000)
    assert np.array_equal(lininterp(xa, np.sin(xa), xva), _fortran_lininterp(xa, np.sin(xa), xva)[0])


def test_read_dsep_drops_truncated_block():
    block = lambda age, n: (["#NPTS MIX-LEN  OVERSH AGE(MYR)", f"#{n} 1.0 0.0 {age}.  0.25 0.004 0.003 -0.76 0.2",
                             "#EEP   MASS    F606W    F606W-F814W"] +
                            [f"{k + 1:3d}  {0.3 + 0.01 * k:.6f} {9.0 - 0.1 * k:.3f}  {1.2 - 0.01 * k:.3f}" for k in range(n)])
    lines = block(6000, 5) + ["", ""] + block(6200, 5)
    with tempfile.NamedTemporaryFile("w", suffix=".iso", delete=False) as fh:
        fh.write("\n".join(lines[:-2]) + "\n 4   0.33")     # second block cut mid-way
    d = read_dsep(fh.name)
    os.remove(fh.name)
    assert sorted(set(d.isochrone_age_yr)) == [6.0e9] and len(d) == 5
    assert np.allclose(d.F814W, d.F606W - np.array([1.2 - 0.01 * k for k in range(5)]))


def test_error_draw():
    e = BandErrors(np.array([100.0, 1e9]), np.array([-14.0, 0.0]), np.array([[1.0, 0.5], [0.0, 1.0]]),
                   {(0, 0): np.array([-0.02, 0.0, 0.01]), (0, 1): np.array([0.1, 0.2]),
                    (1, 0): np.array([0.0]), (1, 1): np.array([-0.05, 0.05])})
    rs = np.random.RandomState(2)
    r, m = np.full(20000, 50.0), np.full(20000, -13.0)        # radial bin 0, magnitude bin 1
    ok, err = e.draw(rs, r, m)
    assert abs(ok.mean() - 0.5) < 0.02 and err[ok].min() >= 0.1 and err[ok].max() <= 0.2
    ok, _ = e.draw(rs, np.full(1000, 500.0), np.full(1000, -15.0))   # bin (1, 0): completeness 0
    assert not ok.any()


def test_reproduces_martins_error_files():
    if not os.path.isdir(MARTIN):
        print("Martin's 47 Tuc files not found; skipped")
        return
    fit = pd.read_csv(f"{MARTIN}/47Tuc_data/47Tuc_fitstars.dat")
    r = fit[(fit.v >= 15.02) & (fit.v <= 19.02)]
    ast = ArtificialStars(f"{MARTIN}/47Tuc_data/47Tucartstars.dat", fit, {"F606W": 30.5742, "F814W": 29.6777},
                          {"F606W": (r.v.min(), r.v.max()), "F814W": (r.i.min(), r.i.max())},
                          {"F606W": 16.48, "F814W": 15.55})
    for b, f in (("F606W", "V"), ("F814W", "I")):
        bd = pd.read_csv(f"{MARTIN}/Photometry/47Tuc_inputfiles/{f}errBoundary.dat", sep=r"\s+", skiprows=3,
                         header=None, names=["rad", "mag", "n", "comp"])
        order = [(i, j) for j in range(12) for i in range(10)]
        assert [len(ast.band[b].errors[k]) for k in order] == bd.n.tolist()
        assert np.allclose([ast.band[b].completeness[k] for k in order], bd.comp.values, atol=1e-6)
        for q, k in enumerate(order):
            ref = np.sort(np.loadtxt(f"{MARTIN}/Photometry/47Tuc_inputfiles/{f}err{q + 1:03d}s.dat", comments="#", ndmin=1))
            assert np.allclose(ref, ast.band[b].errors[k], atol=1e-6)
    d = np.loadtxt(f"{MARTIN}/Photometry/47Tuc_inputfiles/Distance.dat", comments="#")
    assert np.allclose(np.sort(d), ast.dist, atol=1e-4)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(name, "OK")
