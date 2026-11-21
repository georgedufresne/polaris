#!/usr/bin/env python3
"""What the best 47 Tuc sCMDs get wrong: observed vs synthetic colour distributions.

For a grid cell of a variant (default: the best cell of each of D's top 3 runs) the
cell's sCMD is rebuilt from its seed, with each synthetic star's drawn radius kept
(an analysis copy of make_scmd_ast that makes the same random draws in the same
order; the replayed KS is checked against the table). Then, in 0.1-mag F606W
slices of the KS window and in radial zones:

  location  median colour offset from the observed ridge (median colour per slice)
  width     robust sigma (1.4826 MAD) of colour about each sample's own median,
            and the excess sqrt(sig_obs^2 - sig_syn^2)
  tails     fraction redder / bluer than the median by more than 0.03 mag
  KS        where |ECDF_obs - ECDF_syn| (the 2D KS differences) is largest

Fingerprints: PSF term = centred, excess width at all magnitudes, largest relative to
the AST errors in the outskirts; saturation = excess only brighter than ~16.8;
multiple populations = SGB broadening and MS broadening growing fainter; binaries =
red tail; isochrone shape = median offsets drifting along the sequence.

    decompose.py [VARIANT] [RUN:AGE:DM:EBV ...]
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d

F = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(F, "code"))
from isochronetoolbox.scmd.ast import centre                                   # noqa: E402
from isochronetoolbox.scmd.pipeline import A_HST, Setup, apparent, context, select_age  # noqa: E402
from isochronetoolbox.scmd.stats import ks_2d                                   # noqa: E402

rob = lambda q: 1.4826 * np.median(np.abs(q - np.median(q))) if len(q) > 5 else np.nan
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
C_OBS, C_SYN = "#2a78d6", "#eb6834"
ZONES = [(0, 600, "#6da7ec"), (600, 1200, "#256abf"), (1200, 1e9, "#0d366b")]


def scmd_with_radii(rs, mass, phot, ast, n_target, v_window, pad, imf_slope, binary_fraction,
                    martin_sign, rcut, errors):
    """make_scmd_ast (cut="ridge", no extra sigma), also returning each star's radius."""
    V, I = phot["F606W"], phot["F814W"]
    o = np.argsort(mass)
    m, Vs, Is = mass[o], V[o], I[o]
    lo, hi = v_window
    win = np.where((Vs >= lo - pad) & (Vs <= hi + pad))[0]
    a, b = win[0], win[-1]
    m_lo, m_hi = m[a:b + 1].min(), m[a:b + 1].max()
    fV = interp1d(m, Vs, kind="linear", fill_value="extrapolate")
    fI = interp1d(m, Is, kind="linear", fill_value="extrapolate")
    e = imf_slope + 1.0
    sign = 1.0 if martin_sign else -1.0

    def draw(n):
        r = rs.random_sample(n)
        ms = (r * (m_hi ** e - m_lo ** e) + m_lo ** e) ** (1 / e)
        v, i = fV(ms), fI(ms)
        if binary_fraction > 0:
            bb = rs.random_sample(n) <= binary_fraction
            m2 = (0.5 + 0.5 * rs.random_sample(int(bb.sum()))) * ms[bb]
            v[bb] = -2.5 * np.log10(10 ** (-0.4 * v[bb]) + 10 ** (-0.4 * fV(m2)))
            i[bb] = -2.5 * np.log10(10 ** (-0.4 * i[bb]) + 10 ** (-0.4 * fI(m2)))
        rad = ast.radii(rs, n, v)
        if errors in ("nearest", "smooth"):
            ok, ev, ei = ast.nearest.draw(rs, rad, v - ast.zp["F606W"], i - ast.zp["F814W"])
        else:
            okv, ev = ast.band["F606W"].draw(rs, rad, v - ast.zp["F606W"])
            oki, ei = ast.band["F814W"].draw(rs, rad, i - ast.zp["F814W"])
            ok = okv & oki
        return v + sign * ev, i + sign * ei, ok, rad

    def fill(select, n_need):
        out, n_have, drawn, batch = [], 0, 0, max(2 * n_need, 1000)
        while n_have < n_need and drawn < 50.0 * n_need:
            v, i, keep, rad = draw(batch)
            keep &= select(v, i)
            out.append(np.column_stack([v[keep], i[keep], rad[keep]]))
            n_have += int(keep.sum())
            drawn += batch
            rate = max(n_have, 1) / drawn
            batch = int(min(max(1000, 1.2 * (n_need - n_have) / rate), 2e6))
        return np.vstack(out)

    pool = fill(lambda v, i: rcut.in_trace_box(v - i, v), rcut.n_trace)
    t = pool[:rcut.n_trace]
    ridge = rcut.trace(t[:, 0] - t[:, 1], t[:, 0])

    def select(v, i):
        return rcut.select(v - i, v, ridge)
    cat = pool[select(pool[:, 0], pool[:, 1])]
    if len(cat) < n_target:
        cat = np.vstack([cat, fill(select, n_target - len(cat))])
    return cat[:n_target]


def slices(col, mag, rad, ridge_fn, edges):
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (mag >= lo) & (mag < hi)
        c = col[s]
        if len(c) < 30:
            rows.append((np.nan,) * 5)
            continue
        med = np.median(c)
        rows.append((med - ridge_fn((lo + hi) / 2), rob(c), np.mean(c - med > 0.03), np.mean(c - med < -0.03), len(c)))
    return np.array(rows)


def main():
    variant = sys.argv[1] if len(sys.argv) > 1 else "D_zp_mix_nearest"
    s = Setup(**json.load(open(os.path.join(F, variant, "setup.json"))))
    if len(sys.argv) > 2:
        cells = [tuple(float(x) for x in c.split(":")) for c in sys.argv[2:]]
    else:
        per = pd.read_csv(os.path.join(F, "variant_runs.csv"))
        per = per[per.variant == variant].sort_values("ks").head(3)
        cells = [(r.run, r.age, r.dm, r.ebv) for r in per.itertuples()]
    ctx = context(s)
    obs_df = ctx._model_inputs(s, None)                    # ridge-selected observed stars, with x, y
    oc, om = obs_df[s.obs_color].values, obs_df[s.obs_mag].values
    (c0, c1), (m0, m1) = s.ks_window
    inw = (oc >= c0) & (oc <= c1) & (om >= m0) & (om <= m1)
    oc, om = oc[inw], om[inw]
    xc, yc = centre(obs_df.x.values, obs_df.y.values)
    orad = np.hypot(obs_df.x.values - xc, obs_df.y.values - yc)[inw]
    edges = np.round(np.arange(m0, m1 + 1e-9, 0.1), 3)
    mids = 0.5 * (edges[:-1] + edges[1:])
    ridge_med = np.array([np.median(oc[(om >= lo) & (om < hi)]) for lo, hi in zip(edges[:-1], edges[1:])])
    ridge_fn = interp1d(mids, ridge_med, fill_value="extrapolate")
    coarse = np.round(np.arange(m0, m1 + 1e-9, 0.5), 3)
    out = []
    for run, age, dm, ebv in cells:
        run = int(run)
        age_yr, rows = select_age(ctx.isochrones(run), age, s.exact_age)
        coef = {"F606W": 1.0, "F814W": 0.0} if s.reddening == "evi" else A_HST
        phot = apparent(ctx.absolute(rows), dm, ebv, coef, "F606W")
        cat = scmd_with_radii(np.random.RandomState(s.crn_seed), rows["initial_mass"].values, phot, ctx.ast,
                              ctx.n_target, ctx.v_window, s.pad, s.imf_slope, s.binary_fraction,
                              s.martin_sign, ctx.rcut, s.ast_errors)
        sc, sm, srad = cat[:, 0] - cat[:, 1], cat[:, 0], cat[:, 2]
        w = (sc >= c0) & (sc <= c1) & (sm >= m0) & (sm <= m1)
        sc, sm, srad = sc[w], sm[w], srad[w]
        ks, d = ks_2d(ctx.obs3, np.column_stack([sc, sm]), ctx.obs3_cdf, s.ks_percentile)
        tab = pd.read_csv(os.path.join(F, variant, "combined_2d_ks_results.csv"))
        hit = tab[(tab.simulation_id == run) & np.isclose(tab.age_gyr, age) & np.isclose(tab.distance_modulus, dm)
                  & np.isclose(tab.extinction_ebv, ebv)]
        ks_tab = float(hit.ks_statistic.iloc[0]) if len(hit) else np.nan
        print(f"\n=== {variant}  run {run}  {age:.1f} Gyr  mu {dm:.2f}  E {ebv:.3f}:  KS replay {ks:.5f}  table {ks_tab:.5f}"
              f"   n_obs {len(oc)}  n_syn {len(sc)}")
        O = slices(oc, om, orad, ridge_fn, edges)
        S = slices(sc, sm, srad, ridge_fn, edges)
        exc = np.sign(O[:, 1] ** 2 - S[:, 1] ** 2) * np.sqrt(np.abs(O[:, 1] ** 2 - S[:, 1] ** 2))
        print("  F606W   loc_syn-obs  sig_obs  sig_syn  excess   red_obs red_syn  blue_obs blue_syn   [mmag / fractions]")
        for k, mm in enumerate(mids):
            print(f"  {mm:6.2f}  {1000 * (S[k, 0] - O[k, 0]):+8.1f}  {1000 * O[k, 1]:7.1f}  {1000 * S[k, 1]:7.1f}  {1000 * exc[k]:+7.1f}"
                  f"   {O[k, 2]:.3f}  {S[k, 2]:.3f}    {O[k, 3]:.3f}  {S[k, 3]:.3f}")
        print("  radial zones, 0.5-mag slices: sig_obs / sig_syn / excess [mmag]")
        zr = {}
        for lo_r, hi_r, _ in ZONES:
            so, ss = (orad >= lo_r) & (orad < hi_r), (srad >= lo_r) & (srad < hi_r)
            Oz = slices(oc[so], om[so], orad[so], ridge_fn, coarse)
            Sz = slices(sc[ss], sm[ss], srad[ss], ridge_fn, coarse)
            ez = np.sign(Oz[:, 1] ** 2 - Sz[:, 1] ** 2) * np.sqrt(np.abs(Oz[:, 1] ** 2 - Sz[:, 1] ** 2))
            zr[(lo_r, hi_r)] = (Oz, Sz, ez)
            txt = "  ".join(f"{1000 * Oz[k, 1]:5.1f}/{1000 * Sz[k, 1]:5.1f}/{1000 * ez[k]:+5.1f}" for k in range(len(coarse) - 1))
            print(f"   r {lo_r:4.0f}-{min(hi_r, 9999):4.0f}  " + txt)
        print("   (slices " + ", ".join(f"{a:.1f}-{b:.1f}" for a, b in zip(coarse[:-1], coarse[1:])) + ")")
        q95 = np.percentile(d, 95)
        top = d >= q95
        hist, _ = np.histogram(ctx.obs3[top, 1], coarse)
        print("  where |ECDF diff| >= its 95th percentile (the KS value), by F606W: "
              + ", ".join(f"{a:.1f}-{b:.1f}: {h / top.sum():.2f}" for a, b, h in zip(coarse[:-1], coarse[1:], hist)))
        dcol = ctx.obs3[top, 0] - ridge_fn(ctx.obs3[top, 1])
        print(f"  those points' colour offset from the ridge: median {1000 * np.median(dcol):+.1f} mmag,"
              f" quartiles {1000 * np.percentile(dcol, 25):+.1f} .. {1000 * np.percentile(dcol, 75):+.1f}")
        out.append(dict(run=run, age=age, dm=dm, ebv=ebv, ks=ks, O=O, S=S, exc=exc, zr=zr, d=d))

    # figure for the first cell
    c = out[0]
    plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
                         "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                         "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
                         "axes.spines.top": False, "axes.spines.right": False, "font.size": 9})
    fig, ax = plt.subplots(1, 4, figsize=(15, 4.2))
    O, S = c["O"], c["S"]
    ax[0].axhline(0, color=INK2, lw=0.8)
    ax[0].plot(mids, 1000 * (S[:, 0] - O[:, 0]), "-o", color=C_SYN, ms=3, lw=1.5)
    ax[0].set_title("colour offset, synthetic - observed\n(median per slice, mmag)", loc="left", fontsize=9)
    ax[1].plot(mids, 1000 * O[:, 1], "-", color=C_OBS, lw=1.5, label="observed")
    ax[1].plot(mids, 1000 * S[:, 1], "-", color=C_SYN, lw=1.5, label="synthetic")
    ax[1].set_title("colour width (robust sigma, mmag)", loc="left", fontsize=9)
    ax[1].legend(frameon=False)
    for (lo_r, hi_r, colr) in ZONES:
        Oz, Sz, ez = c["zr"][(lo_r, hi_r)]
        cm = 0.5 * (coarse[:-1] + coarse[1:])
        lab = f"r {lo_r:.0f}-{hi_r:.0f} px" if hi_r < 1e8 else f"r > {lo_r:.0f} px"
        ax[2].plot(cm, 1000 * ez, "-o", color=colr, ms=3, lw=1.5, label=lab)
    ax[2].axhline(0, color=INK2, lw=0.8)
    ax[2].axhline(7, color=INK2, lw=0.8, ls=":")
    ax[2].annotate("7 mmag (real-star colour map)", (15.1, 7), xytext=(0, 3), textcoords="offset points", color=INK2, fontsize=7)
    ax[2].set_title("excess colour width, sqrt(obs$^2$ - syn$^2$), mmag\n(negative: model broader)", loc="left", fontsize=9)
    ax[2].legend(frameon=False, fontsize=8)
    ax[3].plot(mids, O[:, 2] - S[:, 2], "-", color="#c0392b" if False else C_SYN, lw=1.5, label="red tail (> +0.03)")
    ax[3].plot(mids, O[:, 3] - S[:, 3], "--", color=C_OBS, lw=1.5, label="blue tail (< -0.03)")
    ax[3].axhline(0, color=INK2, lw=0.8)
    ax[3].set_title("tail fraction, observed - synthetic", loc="left", fontsize=9)
    ax[3].legend(frameon=False, fontsize=8)
    for a_ in ax:
        a_.set_xlabel("F606W")
        for x_, t_ in ((16.5, "short/deep"), (17.02, "SGB"), (17.45, "TO")):
            a_.axvline(x_, color=GRID, lw=1.0, zorder=0)
    fig.suptitle(f"{variant}: run {c['run']}, {c['age']:.1f} Gyr, mu {c['dm']:.2f}, E {c['ebv']:.3f}, KS {c['ks']:.4f}"
                 "   (grey lines: short/deep transition 16.5, SGB 17.0, turnoff 17.45)", x=0.01, ha="left", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    p = os.path.join(F, f"decompose_{variant}_{c['run']}.png")
    fig.savefig(p, dpi=150)
    print("\n" + p)


if __name__ == "__main__":
    main()
