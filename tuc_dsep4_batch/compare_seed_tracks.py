#!/usr/bin/env python3
"""Compare two DSEP4 tracks of the same star run with different numerics
(starting seed or correction damping): log L, log Teff, central H at common ages.

    compare_seed_tracks.py TRACK_A.gz TRACK_B.gz [--from-age 1e8]

Prints the largest differences after the pre-main-sequence transient
(ages >= --from-age) and at a few fixed ages, so a numerics-only change can be
told apart from one that alters the evolution.
"""
import argparse
import gzip

import numpy as np


def read(path):
    lines = gzip.open(path, "rt").read().split("\n")
    h = next(i for i, l in enumerate(lines) if "star_age" in l.split())
    cols = lines[h].split()
    rows = [l.split() for l in lines[h + 1:] if len(l.split()) == len(cols)]
    a = np.array(rows, dtype=float)
    return {c: a[:, i] for i, c in enumerate(cols)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--from-age", type=float, default=1e8)
    o = ap.parse_args()
    A, B = read(o.a), read(o.b)
    lo = max(o.from_age, A["star_age"][0], B["star_age"][0])
    hi = min(A["star_age"][-1], B["star_age"][-1])
    print(f"A: {len(A['star_age'])} models, ages {A['star_age'][0]:.3e}-{A['star_age'][-1]:.3e}")
    print(f"B: {len(B['star_age'])} models, ages {B['star_age'][0]:.3e}-{B['star_age'][-1]:.3e}")
    if hi <= lo:
        print("no common age range past the cut")
        return
    ages = np.geomspace(lo, hi, 400)
    for c in ("log_L", "log_Teff", "center_h1"):
        d = np.interp(ages, A["star_age"], A[c]) - np.interp(ages, B["star_age"], B[c])
        k = np.argmax(np.abs(d))
        print(f"  max |A-B| {c:10s} = {abs(d[k]):.2e}  (at {ages[k]:.3e} yr)   median |A-B| = {np.median(np.abs(d)):.2e}")
    for t in (1e8, 1e9, 5e9, 1e10, 1.6e10):
        if lo <= t <= hi:
            vals = [f"{c}: {np.interp(t, A['star_age'], A[c]):.5f} vs {np.interp(t, B['star_age'], B[c]):.5f}"
                    for c in ("log_L", "log_Teff", "center_h1")]
            print(f"  {t:.1e} yr  " + "   ".join(vals))


if __name__ == "__main__":
    main()
