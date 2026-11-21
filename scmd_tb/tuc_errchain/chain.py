#!/usr/bin/env python3
"""From Martin's error model to ours, one change at a time, on 47 Tuc (our ridge cut and window throughout).

  python chain.py fit  VARIANT RUN   -> fits/VARIANT/RUN.csv   KS on a local (age, mu, E) grid, CRN draws
  python chain.py null VARIANT       -> nulls/VARIANT.csv      bootstrap null (fresh references, resampled
                                                               ASTs) at the variant's best cell
All variants use the literature PDMF slope (-0.45) and binary fraction (0.009, q > 0.5)."""
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd

E = "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_errchain"
sys.path.insert(0, f"{E}/code")
from isochronetoolbox.scmd import Setup, context                   # noqa: E402
from isochronetoolbox.scmd.pipeline import _window                  # noqa: E402
from isochronetoolbox.scmd.stats import ks_2d                       # noqa: E402
from isochronetoolbox.scmd import bootstrap as BS                   # noqa: E402

LIT = {"imf_slope": -0.45, "binary_fraction": 0.009}
ZPK = float(__import__("os").environ.get("TUC_ZPK", "0.0"))
RAW = {"ast_errors": "nearest", "radial": "corrected", "radial_source": "raw"}
VARIANTS = {
    "V0_martin":    {"ast_errors": "independent", "radial": "observed"},        # MakeFakeCMD errors + radii
    "V1_joint":     {"ast_errors": "nearest", "radial": "observed"},            # + joint errors & completeness
    "V2_weighted":  {"ast_errors": "nearest", "radial": "corrected"},           # + 1/c radii (production)
    "V3_rawradial": dict(RAW),                                                  # + radii from before the cut
    "V4_floor":     dict(RAW, radial_cfloor=0.02),
    "S_k15":        dict(RAW, ast_k=15),
    "S_k60":        dict(RAW, ast_k=60),
    "S_scale_half": dict(RAW, ast_scale_mag=0.05, ast_scale_lnr=0.045),
    "S_scale_dbl":  dict(RAW, ast_scale_mag=0.2, ast_scale_lnr=0.18),
    "S_split":      dict(RAW, ast_split=True),
    "S_mag_half":   dict(RAW, ast_scale_mag=0.05),            # neighbour sets depend only on the mag/radius
    "S_mag_dbl":    dict(RAW, ast_scale_mag=0.2),             # scale ratio: change the magnitude scale alone
    "V5_zpcorr":    dict(RAW, zp_correct=True, zp_ref=(17.9, 19.0)),
    "V6_zpsplit":   dict(RAW, zp_correct=True, zp_ref=(17.9, 19.0), zp_mag_frac=ZPK),
}
AGES = np.round(np.arange(11.0, 13.6001, 0.1), 2)
MUS = np.round(np.arange(13.30, 13.4801, 0.035), 3)
ES = np.round(np.arange(0.000, 0.0501, 0.005), 3)

mode, VAR = sys.argv[1], sys.argv[2]
d = {k: v for k, v in json.load(open(f"{E}/setup.json")).items() if k in Setup.__dataclass_fields__}
d.update(LIT)
d.update(VARIANTS[VAR])
S = Setup(**d)
CTX = context(S)
MAKE = {}


def one(job):
    run, age, mu, e = job
    if (run, age) not in MAKE:
        MAKE[(run, age)] = CTX.model(run, age)[1]
    p = CTX.points(MAKE[(run, age)](np.random.RandomState(S.crn_seed), mu, e))
    p = p[_window(p, S.ks_window)]
    return run, age, mu, e, ks_2d(CTX.obs3, p, CTX.obs3_cdf)[0]


def null_slice(args):
    cell, ids = args
    rows = BS.realisations(CTX, cell["run"], cell["age"], cell["mu"], cell["evi"], ids, fresh_reference=True)
    return [r["ks_statistic"] for r in rows]


if __name__ == "__main__":
    if mode == "fit":
        run = int(sys.argv[3])
        os.makedirs(f"{E}/fits/{VAR}", exist_ok=True)
        jobs = [(run, a, m, e) for a in AGES for m in MUS for e in ES]
        with Pool(8) as pool:
            res = pool.map(one, jobs, chunksize=6)
        df = pd.DataFrame(res, columns=["run", "age", "mu", "evi", "ks"])
        df.to_csv(f"{E}/fits/{VAR}/{run}.csv", index=False)
        b = df.loc[df.ks.idxmin()]
        print(f"{VAR} run {run}: min KS {b.ks:.5f} at {b.age} Gyr, mu {b.mu}, E {b.evi}")
    elif mode == "null":
        n = int(sys.argv[3]) if len(sys.argv) > 3 else 240
        f = pd.concat([pd.read_csv(f"{E}/fits/{VAR}/{x}") for x in os.listdir(f"{E}/fits/{VAR}")], ignore_index=True)
        b = f.loc[f.ks.idxmin()]
        cell = {"run": int(b.run), "age": float(b.age), "mu": float(b.mu), "evi": float(b.evi)}
        with Pool(8) as pool:
            ks = np.concatenate(pool.map(null_slice, [(cell, list(range(k, n, 8))) for k in range(8)]))
        os.makedirs(f"{E}/nulls", exist_ok=True)
        pd.DataFrame({"ks": ks}).assign(**cell).to_csv(f"{E}/nulls/{VAR}.csv", index=False)
        print(f"{VAR} null at {cell}: n {len(ks)} mean {ks.mean():.5f} std {ks.std():.5f} "
              f"3sigma {ks.mean() + 3 * ks.std():.5f}; best KS {b.ks:.5f} = {(b.ks - ks.mean()) / ks.std():.1f} sigma")
