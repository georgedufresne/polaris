#!/usr/bin/env python3
"""KS (95th) of one cell under the new model for a given PDMF slope (common random numbers)."""
import itertools, json, sys
from dataclasses import replace
import numpy as np, pandas as pd
O, task = sys.argv[1], int(sys.argv[2])
sys.path.insert(0, f"{O}/code")
from isochronetoolbox.scmd import Setup, context
alphas = [-0.84, -0.45]
ages = [11.2, 11.4, 11.6, 11.8]
mus = [13.34, 13.36, 13.38, 13.40, 13.42, 13.44]
evis = [0.005, 0.015, 0.025]
cells = list(itertools.product(alphas, ages))           # one (alpha, age) per task, all (mu, E) inside
alpha, age = cells[task]
S = Setup(**{k: v for k, v in json.load(open(f"{O}/setup.json")).items() if k in Setup.__dataclass_fields__})
ctx = context(replace(S, imf_slope=alpha))
_, make = ctx.model(10099, age)
rows = []
for mu in mus:
    for e in evis:
        cat = make(np.random.RandomState(S.crn_seed), mu, e)
        r = ctx.stage3(cat)
        p = ctx.points(cat)
        rows.append({"alpha": alpha, "age": age, "mu": mu, "evi": e, "ks": r["ks_statistic"],
                     "frac_above_17": float(np.mean(p[(p[:, 1] > 15.02) & (p[:, 1] < 19.02), 1] < 17.0))})
        print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(f"{O}/parts/part_{task:02d}.csv", index=False)
