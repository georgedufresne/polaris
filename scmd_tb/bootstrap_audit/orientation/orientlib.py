"""2D-KS statistics in the four ECDF orientations (audit helper, not part of the toolbox).

Orientation = which quadrant the ECDF counts.  Coordinates are (colour, F606W);
rank2d counts points strictly below-left, so multiplying the columns by the
signs below before ranking/interpolating picks the quadrant:

  bb (+c, +m): bluer  & brighter  -- the pipeline (stats.ks_2d)
  bf (+c, -m): bluer  & fainter   -- M. Ying (vorbin_chi2.KS_2d negates V)
  rb (-c, +m): redder & brighter
  rf (-c, -m): redder & fainter
"""
import numpy as np
from scipy.interpolate import LinearNDInterpolator

ORIENT = [("bb", (1.0, 1.0)), ("bf", (1.0, -1.0)), ("rb", (-1.0, 1.0)), ("rf", (-1.0, -1.0))]
KEYS = [k for k, _ in ORIENT]


def ref_ecdfs(ref, rank2d):
    """Own ECDF of the reference set in each orientation (rank2d / N on transformed coordinates)."""
    return {k: rank2d(ref * np.array(s)) / len(ref) for k, s in ORIENT}


def pct(a, q):
    return float(np.percentile(a, q)) if len(a) else float("nan")


def ying_stat(ref, sim, fi_bf, F_bf, q):
    """vorbin_chi2.KS_2d.evaluate, with the sCMD ECDF already interpolated (fi_bf, orientation bf)."""
    obs_t = ref * np.array([1.0, -1.0])
    fit_t = sim * np.array([1.0, -1.0])
    max_v, min_v = obs_t[:, 1].max(), obs_t[:, 1].min()
    lower = fit_t[fit_t[:, 1] < min_v]
    upper_n = int((fit_t[:, 1] < max_v).sum())
    lower_n = len(lower)
    lower_vi = np.sort(lower[:, 0])
    pred = (fi_bf * len(fit_t) - np.searchsorted(lower_vi, obs_t[:, 0], side="left")) / (upper_n - lower_n)
    mask = (F_bf != 0.0) & (pred != 0.0)
    delta = np.abs(pred[mask] - F_bf[mask])
    k = int(len(delta) * q / 100.0)
    return float(np.partition(delta, k)[k]), int(mask.sum()), upper_n, lower_n


def score(ref, F_ref, sim, rank2d, validate=False):
    """All statistics of one sCMD ``sim`` against ``ref`` (both (colour, F606W), in the KS window).

    Interpolated ECDF: the pipeline's recipe (rank2d/M at the sCMD's points, LinearNDInterpolator
    with fill 0) on coordinates transformed per orientation, identically for both sets.
    Exact ECDF: count of sCMD stars in the quadrant of each reference star (one masked rank2d
    in the bb orientation plus 1D counts; ``validate`` checks this against a literal masked
    rank2d in every orientation).
    """
    N, M = len(ref), len(sim)
    out = {"n_ref": N, "n_sim": M}
    arr = {}
    both = np.vstack([sim, ref])
    mask = np.r_[np.ones(M, bool), np.zeros(N, bool)]
    n_bb = rank2d(both, mask)[M:].astype(float)
    sc, sm = np.sort(sim[:, 0]), np.sort(sim[:, 1])
    n_c = np.searchsorted(sc, ref[:, 0], side="left").astype(float)     # sCMD stars strictly bluer
    n_m = np.searchsorted(sm, ref[:, 1], side="left").astype(float)     # strictly brighter
    out["cross_ties_c"] = int((np.searchsorted(sc, ref[:, 0], side="right") - n_c).sum())
    out["cross_ties_m"] = int((np.searchsorted(sm, ref[:, 1], side="right") - n_m).sum())
    n_ex = {"bb": n_bb, "bf": n_c - n_bb, "rb": n_m - n_bb, "rf": M - n_c - n_m + n_bb}
    # 1D marginal CDFs (strict) at the reference points
    F1c = np.searchsorted(np.sort(ref[:, 0]), ref[:, 0], side="left") / N
    F1m = np.searchsorted(np.sort(ref[:, 1]), ref[:, 1], side="left") / N
    d1c, d1m = F1c - n_c / M, F1m - n_m / M
    out.update(c1d_p95=pct(np.abs(d1c), 95), c1d_p99=pct(np.abs(d1c), 99), c1d_max=float(np.abs(d1c).max()),
               m1d_p95=pct(np.abs(d1m), 95), m1d_p99=pct(np.abs(d1m), 99), m1d_max=float(np.abs(d1m).max()))
    arr["d1c"], arr["d1m"] = d1c, d1m
    for k, s in ORIENT:
        sg = np.array(s)
        st, rt = sim * sg, ref * sg
        fs = rank2d(st) / M
        fi = LinearNDInterpolator(st, fs, fill_value=0)(rt)
        d = F_ref[k] - fi
        dx = F_ref[k] - n_ex[k] / M
        ad, adx = np.abs(d), np.abs(dx)
        z = (F_ref[k] != 0) & (fi != 0)
        out.update({f"{k}_p95": pct(ad, 95), f"{k}_p99": pct(ad, 99), f"{k}_max": float(ad.max()),
                    f"{k}_mean": float(ad.mean()),
                    f"{k}_z_p95": pct(ad[z], 95), f"{k}_z_p99": pct(ad[z], 99),
                    f"{k}_nzero": int((~z).sum()), f"{k}_nfi0": int((fi == 0).sum()),
                    f"{k}_ex_p95": pct(adx, 95), f"{k}_ex_p99": pct(adx, 99), f"{k}_ex_max": float(adx.max())})
        arr[f"d_{k}"], arr[f"dx_{k}"], arr[f"F_{k}"], arr[f"fi_{k}"] = d, dx, F_ref[k], fi
        if validate:
            lit = rank2d(np.vstack([st, rt]), mask)[M:].astype(float)
            out[f"{k}_ex_identity_maxdiff"] = float(np.abs(lit - n_ex[k]).max())
        if k == "bf":
            for q in (95, 99):
                v, nm, un, ln = ying_stat(ref, sim, fi, F_ref["bf"], q)
                out[f"ying_p{q}"] = v
            out["ying_nmask"], out["ying_upper_n"], out["ying_lower_n"] = nm, un, ln
    # Fasano-Franceschini style, two ways: per-point max over the 4 quadrants, then the percentile ...
    for tag, pre in (("", "d_"), ("ex_", "dx_")):
        mx = np.max(np.abs(np.vstack([arr[pre + k] for k in KEYS])), axis=0)
        out[f"ffpt_{tag}p95"], out[f"ffpt_{tag}p99"], out[f"ffpt_{tag}max"] = pct(mx, 95), pct(mx, 99), float(mx.max())
    # ... and the max over the 4 per-orientation statistics
    for st in ("p95", "p99", "max", "ex_p95", "ex_p99", "ex_max", "z_p95", "z_p99"):
        out[f"ffmax_{st}"] = max(out[f"{k}_{st}"] for k in KEYS)
    return out, arr
