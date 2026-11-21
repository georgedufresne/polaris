#!/usr/bin/env python3
"""Stage-3 KS for every (mu, reddening) cell of one (run, age), not only stage 2's best four."""
import json
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

G = sys.argv[1]
RUN, AGE = int(sys.argv[2]), float(sys.argv[3])
sys.path.insert(0, f"{G}/code")
from isochronetoolbox.scmd import Setup, context                   # noqa: E402
from isochronetoolbox.scmd.pipeline import _window                  # noqa: E402
from isochronetoolbox.scmd.stats import ks_2d                       # noqa: E402

S = Setup(**{k: v for k, v in json.load(open(f"{G}/setup.json")).items() if k in Setup.__dataclass_fields__})
CTX = context(S)


def score(item):
    dm, e, sim = item
    js = CTX.stage2(None if sim is None else sim_cat[(dm, e)])
    return dm, e


def ks_one(args):
    dm, e, pts = args
    return dm, e, ks_2d(CTX.obs3, pts, CTX.obs3_cdf)[0]


if __name__ == "__main__":
    age_yr, cells = CTX.grid(RUN, AGE)
    items, js = [], {}
    for dm, e, cat in cells:                      # in the grid's order, from the grid's seed
        if cat is None:
            continue
        r = CTX.stage2(cat)
        js[(round(dm, 3), round(e, 4))] = r[0] if r else np.nan
        p = CTX.points(cat)
        items.append((round(dm, 3), round(e, 4), p[_window(p, S.ks_window)]))
    with Pool(8) as pool:
        res = pool.map(ks_one, items)
    df = pd.DataFrame(res, columns=["mu", "evi", "ks"])
    df["js"] = [js[(m, e)] for m, e in zip(df.mu, df.evi)]
    df["js_rank"] = df.js.rank(method="first").astype(int)
    df.to_csv(f"{G}/all_cells_{RUN}_{AGE}.csv", index=False)
    grid = pd.read_csv(f"{G}/combined_2d_ks_results.csv")
    g = grid[(grid.simulation_id == RUN) & (np.isclose(grid.age_gyr, AGE))]
    print(f"run {RUN}, {AGE} Gyr: {len(df)} cells")
    print(f"  grid (stage-2 best 4): min KS {g.ks_statistic.min():.5f} at mu {g.loc[g.ks_statistic.idxmin(), 'distance_modulus']}, E {g.loc[g.ks_statistic.idxmin(), 'extinction_ebv']}")
    b = df.loc[df.ks.idxmin()]
    print(f"  all 341 cells        : min KS {b.ks:.5f} at mu {b.mu}, E {b.evi} (its stage-2 JS rank {int(b.js_rank)} of {len(df)})")
    print("  10 lowest-KS cells:")
    print(df.nsmallest(10, "ks").to_string(index=False, float_format=lambda v: f"{v:.5f}"))
    print(f"  Spearman(JS, KS) over all cells: {df.js.corr(df.ks, method='spearman'):.3f}")
