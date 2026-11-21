#!/usr/bin/env python3
"""Diagnostic: KS of one run on an (age, mu, E(V-I)) grid that extends E(V-I) below the 0.08 prior floor.
Same model as the grid (m15_grid_d4/setup.json), seeds 7000 + k (not the grid's CRN), output one CSV per run."""
import json, sys
from multiprocessing import Pool
import numpy as np, pandas as pd

G, RUN = sys.argv[1], int(sys.argv[2])
AGES = np.round(np.arange(11.5, 15.0001, 0.5), 2)
MUS = np.round(np.arange(15.42, 15.5801, 0.02), 3)
ES = np.round(np.arange(0.00, 0.1401, 0.02), 3)
sys.path.insert(0, f"{G}/code")
from isochronetoolbox.scmd import Setup, context                   # noqa: E402
from isochronetoolbox.scmd.pipeline import _window                  # noqa: E402
from isochronetoolbox.scmd.stats import ks_2d                       # noqa: E402
S = Setup(**{k: v for k, v in json.load(open(f"{G}/setup.json")).items() if k in Setup.__dataclass_fields__})
CTX = context(S)
MAKE = {}


def one(args):
    k, age, mu, e = args
    if age not in MAKE:
        MAKE[age] = CTX.model(RUN, age)[1]
    p = CTX.points(MAKE[age](np.random.RandomState(7000 + k), mu, e))
    p = p[_window(p, S.ks_window)]
    return RUN, age, mu, e, ks_2d(CTX.obs3, p, CTX.obs3_cdf)[0]


if __name__ == "__main__":
    jobs = [(k, a, m, e) for k, (a, m, e) in enumerate((a, m, e) for a in AGES for m in MUS for e in ES)]
    with Pool(8) as pool:
        res = pool.map(one, jobs, chunksize=4)
    df = pd.DataFrame(res, columns=["run", "age", "mu", "evi", "ks"])
    df.to_csv(f"{G}/escan/escan_{RUN}.csv", index=False)
    b = df.loc[df.ks.idxmin()]
    print(f"run {RUN}: min KS {b.ks:.5f} at {b.age} Gyr, mu {b.mu}, E {b.evi}")
