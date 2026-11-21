"""Tests for isochronetoolbox.scmd.ridge.

  - the vectorised DP equals the original loop (transcribed below) to the bit
  - the vectorised snap equals the original scikit-learn loop
  - the binned KDE gives the exact KDE's ridge
  - RidgeCut keeps the stars near a known line and fixes its scales from the data

Run directly (``python tests/scmd/test_ridge.py``) or under pytest.
"""

from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(os.path.dirname(HERE)))

from isochronetoolbox.scmd import ridge as R  # noqa: E402


def _cmd(n=6000, seed=1):
    """A turnoff-shaped sequence: main sequence, subgiant branch, giant branch."""
    rs = np.random.RandomState(seed)
    t = rs.uniform(0, 1, n)
    mag = 19.0 - 4.0 * t
    color = 0.52 + 0.4 * (t - 0.4) ** 2 + 0.12 / (1 + np.exp(-(t - 0.55) / 0.02))
    return color + rs.normal(0, 0.012, n), mag + rs.normal(0, 0.01, n)


def _dp_loop(cost, cols, mags_grid, dc_max, lam_slope, lam_curv, bias):
    """ridge_via_dp's dynamic programme as originally written (per-cell loop)."""
    n_mag, n_col = cost.shape
    J = np.full((n_mag, n_col), np.inf)
    prev = np.full((n_mag, n_col), -1, dtype=int)
    J[0, :] = cost[0, :]
    max_jump = max(1, int(np.ceil(dc_max / (cols[1] - cols[0]))))
    for i in range(1, n_mag):
        for j in range(n_col):
            jlo, jhi = max(0, j - max_jump), min(n_col, j + max_jump + 1)
            src = np.arange(jlo, jhi)
            slope = (cols[j] - cols[src]) / (mags_grid[i - 1] - mags_grid[i] + 1e-10)
            cand = J[i - 1, jlo:jhi] + cost[i, j] + bias[i, j] + lam_slope * slope ** 2
            if i >= 2:
                jp = prev[i - 1, jlo:jhi]
                ok = (jp >= 0) & (jp < n_col)
                if ok.any():
                    curv = np.zeros(len(src))
                    curv[ok] = cols[j] - 2 * cols[src][ok] + cols[jp[ok]]
                    cand = cand + lam_curv * curv ** 2
            kb = int(np.argmin(cand))
            J[i, j] = cand[kb]
            prev[i, j] = src[kb]
    j = int(np.argmin(J[-1, :]))
    path = [j]
    for i in range(n_mag - 1, 0, -1):
        j = prev[i, j]
        if j < 0 or j >= n_col:
            break
        path.append(j)
    path = path[::-1]
    return np.column_stack([cols[np.array(path)], mags_grid[:len(path)]])


def test_dp_matches_loop():
    c, m = _cmd()
    bw = R.bandwidth(c, m)
    for red in (None, 16.5):
        n_mag, n_col = 60, 120
        got = R.ridge_via_dp(c, m, n_mag=n_mag, n_col=n_col, dc_max=0.035, lam_slope=3.0, lam_curv=15.0,
                             red_bias_after=red, lam_red=1.5, bw=bw, kde="binned")
        cols = np.linspace(c.min() - 0.05, c.max() + 0.05, n_col)
        mags_grid = np.linspace(m.max() + 0.05, m.min() - 0.05, n_mag)
        cost = -R._log_density(c, m, cols, mags_grid, bw, "binned")
        bias = np.zeros_like(cost)
        if red is not None:
            j0 = np.where(mags_grid <= red)[0][0]
            bias[j0:, :] -= (1.5 * np.maximum(0, cols - cols.mean()))[np.newaxis, :]
        ref = _dp_loop(cost, cols, mags_grid, 0.035, 3.0, 15.0, bias)
        assert np.array_equal(got, ref)


def test_snap_matches_sklearn():
    from sklearn.decomposition import PCA
    from sklearn.neighbors import NearestNeighbors
    c, m = _cmd()
    ridge = R.ridge_via_dp(c, m, n_mag=60, n_col=120, bw=R.bandwidth(c, m), kde="binned")
    X = np.column_stack([c, m])
    nbrs = NearestNeighbors(n_neighbors=80).fit(X)
    ref = ridge.copy()
    for _ in range(2):
        for i, point in enumerate(ref):
            dist, idx = nbrs.kneighbors([point], n_neighbors=80)
            nb = X[idx[0]]
            normal = PCA(n_components=2).fit(nb - nb.mean(axis=0)).components_[1]
            proj = np.dot(nb - point, normal)
            w = np.exp(-dist[0] / np.mean(dist[0]))
            o = np.argsort(proj)
            cs = np.cumsum(w[o])
            k = min(np.searchsorted(cs, cs[-1] * 0.5), len(proj) - 1)
            ref[i] = point + proj[o[k]] * normal
    got = R.snap_to_local_ridge(ridge, c, m, k_neighbors=80, n_iter=2)
    assert np.allclose(got, ref, rtol=0, atol=1e-10)


def test_binned_kde_gives_exact_ridge():
    c, m = _cmd()
    bw = R.bandwidth(c, m)
    ex, bi = R.Tracer(kde="exact")(c, m, bw), R.Tracer(kde="binned")(c, m, bw)
    sc = (np.ptp(c), np.ptp(m))
    assert np.max(R.Ridge(ex[:, 0], ex[:, 1], sc).distance(bi[:, 0], bi[:, 1])) < 0.005


def test_ridge_cut():
    c, m = _cmd(20000)
    rs = np.random.RandomState(2)
    c[:500] += rs.uniform(0.1, 0.3, 500)                     # a red tail the cut must remove
    rc, keep = R.RidgeCut.from_observed(R.Tracer(), c, m, (15.2, 18.8), (0.3, 1.2), 0.05)
    box, tb = rc.in_box(c, m), rc.in_trace_box(c, m)
    assert rc.n_trace == tb.sum() > box.sum() and np.isclose(rc.scales[0], np.ptp(c[box]))
    assert not keep[~box].any()
    assert keep[:500].mean() < 0.1 and keep[500:][box[500:]].mean() > 0.97
    again = rc.trace(c[tb], m[tb])                          # the same stars give the same ridge
    assert np.array_equal(again.distance(c, m), rc.ridge.distance(c, m))
    ex = ((0.6, 0.7), (15.5, 16.5))                         # an excluded box: out of the trace and the cut
    rc2, keep2 = R.RidgeCut.from_observed(R.Tracer(), c, m, (15.2, 18.8), (0.3, 1.2), 0.05, exclude=[ex])
    inside = (c > 0.6) & (c < 0.7) & (m > 15.5) & (m < 16.5)
    assert inside.any() and not keep2[inside].any() and rc2.n_trace == rc.n_trace - np.sum(inside & tb)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(name, "OK")
