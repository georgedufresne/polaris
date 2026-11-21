#!/usr/bin/env python3
"""Data vs fresh realisations of one cell (stage-3 KS), under the setup in OUT."""
import json, sys
import numpy as np, pandas as pd
OUT, cell, n = sys.argv[1], sys.argv[2], int(sys.argv[3])
sys.path.insert(0, f"{OUT}/code")
from isochronetoolbox.scmd import Setup, context
S = Setup(**{k: v for k, v in json.load(open(f"{OUT}/setup.json")).items() if k in Setup.__dataclass_fields__})
ctx = context(S)
r, a, dm, e = cell.split(":"); r, a, dm, e = int(r), float(a), float(dm), float(e)
_, make = ctx.model(r, a)
rows = []
for i in range(n):
    res = ctx.stage3(make(np.random.RandomState(500000 + i), dm, e)); res["seed"] = 500000 + i; rows.append(res)
    print(i, res["ks_statistic"], flush=True)
pd.DataFrame(rows).to_csv(f"{OUT}/data_vs_model.csv", index=False)
