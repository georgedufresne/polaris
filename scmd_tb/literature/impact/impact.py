#!/usr/bin/env python3
"""KS at fixed best-fit cells with the old vs literature PDMF slope / binary fraction, several seeds each."""
import json, sys
from multiprocessing import Pool
import numpy as np, pandas as pd

G = sys.argv[1]; CELLS = json.loads(sys.argv[2]); VARIANTS = json.loads(sys.argv[3]); NSEED = 6
sys.path.insert(0, f"{G}/code")
from isochronetoolbox.scmd import Setup, context                   # noqa: E402
from isochronetoolbox.scmd.pipeline import _window                  # noqa: E402
from isochronetoolbox.scmd.stats import ks_2d                       # noqa: E402
S = Setup(**{k: v for k, v in json.load(open(f"{G}/setup.json")).items() if k in Setup.__dataclass_fields__})
CTX = context(S)


def one(job):
    (run, age, mu, e), (alpha, fb), seed = job
    object.__setattr__(CTX.s, "imf_slope", alpha); object.__setattr__(CTX.s, "binary_fraction", fb)
    make = CTX.model(run, age)[1]
    p = CTX.points(make(np.random.RandomState(31000 + seed), mu, e))
    return run, age, mu, e, alpha, fb, seed, ks_2d(CTX.obs3, p[_window(p, S.ks_window)], CTX.obs3_cdf)[0]


if __name__ == "__main__":
    jobs = [(tuple(c), tuple(v), s) for c in CELLS for v in VARIANTS for s in range(NSEED)]
    with Pool(8) as pool:
        r = pd.DataFrame(pool.map(one, jobs), columns=["run", "age", "mu", "evi", "alpha", "fbin", "seed", "ks"])
    r.to_csv(f"{sys.argv[4]}", index=False)
    t = r.groupby(["run", "age", "mu", "evi", "alpha", "fbin"]).ks.agg(["mean", "std"]).reset_index()
    print(t.to_string(index=False, float_format=lambda v: f"{v:.5f}"))
