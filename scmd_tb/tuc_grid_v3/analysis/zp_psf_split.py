#!/usr/bin/env python3
"""KS of the v2 and v3 best cells under the four combinations of zero-point
correction (data) and PSF scatter (model), grid's frozen code."""
import json, sys
from multiprocessing import Pool
import numpy as np, pandas as pd
B = "/dartfs-hpc/scratch/f007058/scmd_tb"
sys.path.insert(0, f"{B}/tuc_grid_v2/code")
CELLS = [(10005, 11.1, 13.44, 0.030), (10035, 11.4, 13.38, 0.000), (10030, 12.3, 13.34, 0.005),
         (10091, 11.4, 13.44, 0.035), (10070, 11.1, 13.43, 0.030), (10046, 12.5, 13.32, 0.000)]
PSF = [[16.0, 0.0074], [17.8, 0.0031], [18.2, 0.004], [18.6, 0.0049], [18.9, 0.0059]]


def one(combo):
    zp, psf = combo
    from isochronetoolbox.scmd import Setup, context, regenerate
    saved = json.load(open(f"{B}/tuc_grid_v2/setup.json"))
    saved.update(zp_correct=zp, extra_sigma_table=PSF if psf else [])
    ctx = context(Setup(**{k: v for k, v in saved.items() if k in Setup.__dataclass_fields__}))
    return combo, [ctx.stage3(regenerate(ctx, *c))["ks_statistic"] for c in CELLS]


with Pool(4) as p:
    res = dict(p.map(one, [(True, True), (False, True), (True, False), (False, False)]))
null = {True: pd.read_csv(f"{B}/tuc_grid_v2/bootstrap_2d_ks_results.csv").ks_statistic.values,
        False: pd.read_csv(f"{B}/tuc_grid_v3/bootstrap_2d_ks_results.csv").ks_statistic.values}
rows = []
for k, c in enumerate(CELLS):
    r = {"cell": f"{c[0]} {c[1]:.1f} Gyr mu {c[2]:.2f} E {c[3]:.3f}", "from": "v2 top3" if k < 3 else "v3 top3"}
    for (zp, psf), ks in res.items():
        n = null[psf]
        r[f"zp={'on' if zp else 'off'} psf={'on' if psf else 'off'}"] = f"{ks[k]:.4f} (p {(n >= ks[k]).mean():.3f})"
    rows.append(r)
df = pd.DataFrame(rows); print(df.to_string(index=False))
for psf in (True, False):
    n = null[psf]; print(f"null psf={'on' if psf else 'off'}: {n.mean():.4f} +- {n.std():.4f}, 3 sigma {n.mean() + 3 * n.std():.4f}")
df.to_csv(f"{B}/tuc_grid_v3/analysis/zp_psf_split.csv", index=False)
