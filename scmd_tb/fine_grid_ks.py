#!/usr/bin/env python3
"""KS of one run/age on a fine (mu, E(V-I)) grid: is the best colour offset between the grid's 0.005 steps?"""
import json
import sys
from multiprocessing import Pool

import numpy as np

G = sys.argv[1]
RUN, AGE = int(sys.argv[2]), float(sys.argv[3])
MUS = np.round(np.arange(float(sys.argv[4]), float(sys.argv[5]) + 1e-9, float(sys.argv[6])), 4)
ES = np.round(np.arange(float(sys.argv[7]), float(sys.argv[8]) + 1e-9, float(sys.argv[9])), 4)
import os
sys.path.insert(0, f"{G}/code_boot" if os.path.isdir(f"{G}/code_boot") else f"{G}/code")
from isochronetoolbox.scmd import Setup, context                   # noqa: E402
from isochronetoolbox.scmd.pipeline import _window                  # noqa: E402
from isochronetoolbox.scmd.stats import ks_2d                       # noqa: E402

S = Setup(**{k: v for k, v in json.load(open(f"{G}/setup.json")).items() if k in Setup.__dataclass_fields__})
CTX = context(S)
_, MAKE = CTX.model(RUN, AGE)


def one(args):
    k, mu, e = args
    p = CTX.points(MAKE(np.random.RandomState(5000 + k), mu, e))
    p = p[_window(p, S.ks_window)]
    return mu, e, ks_2d(CTX.obs3, p, CTX.obs3_cdf)[0]


if __name__ == "__main__":
    jobs = [(k, mu, e) for k, (mu, e) in enumerate((m, x) for m in MUS for x in ES)]
    with Pool(8) as pool:
        res = pool.map(one, jobs)
    tab = {(m, e): v for m, e, v in res}
    print(f"run {RUN}, {AGE} Gyr: KS on the fine grid (rows mu, columns E(V-I))")
    print("  mu \\ E " + " ".join(f"{e:7.3f}" for e in ES))
    for m in MUS:
        print(f"  {m:7.3f} " + " ".join(f"{tab[(m, e)]:7.4f}" for e in ES))
    b = min(res, key=lambda r: r[2])
    print(f"minimum {b[2]:.5f} at mu {b[0]}, E(V-I) {b[1]}")
