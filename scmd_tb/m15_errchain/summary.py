#!/usr/bin/env python3
"""Collate the error-model chain: best KS, best (age, mu, E), per-run best age, null, per variant."""
import glob
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

E = "/dartfs-hpc/scratch/f007058/scmd_tb/m15_errchain"
ORDER = ["V0_martin", "V1_joint", "V2_weighted", "V3_rawradial", "V4_floor", "V5_zpcorr", "V6_zpsplit",
         "S_k15", "S_k60", "S_mag_half", "S_mag_dbl", "S_split"]
LABEL = {"V0_martin": "Martin: per-band errors, unweighted radii", "V1_joint": "+ joint errors & completeness",
         "V2_weighted": "+ 1/c radii (production)", "V3_rawradial": "+ radii from before the cut",
         "V4_floor": "  c floor 0.02", "V5_zpcorr": "+ zero-point map (colour only)",
         "V6_zpsplit": "+ zero-point map (both bands)",
         "S_k15": "  k = 15", "S_k60": "  k = 60", "S_mag_half": "  mag scale x 0.5 (radius fixed)",
         "S_mag_dbl": "  mag scale x 2 (radius fixed)", "S_split": "  short/long split"}
rows, prof, perrun = [], {}, {}
for v in ORDER:
    fs = glob.glob(f"{E}/fits/{v}/*.csv")
    if not fs:
        continue
    f = pd.concat([pd.read_csv(p) for p in fs], ignore_index=True)
    b = f.loc[f.ks.idxmin()]
    b8 = f[f.evi >= 0.08 - 1e-9].loc[lambda d: d.ks.idxmin()]
    pr = f.groupby("age").ks.min()
    prof[v] = pr
    pb = f.loc[f.groupby("run").ks.idxmin()].set_index("run")
    perrun[v] = pb
    r = {"variant": v, "runs": f.run.nunique(), "best_ks": b.ks, "run": int(b.run), "age": b.age, "mu": b.mu,
         "evi": b.evi, "best_ks_E>=0.08": b8.ks, "age_E>=0.08": b8.age,
         "median_run_best_age": pb.age.median(), "median_run_best_E": pb.evi.median()}
    nf = f"{E}/nulls/{v}.csv"
    if os.path.exists(nf):
        n = pd.read_csv(nf)
        r.update(null_mean=n.ks.mean(), null_sd=n.ks.std(), null_3sig=n.ks.mean() + 3 * n.ks.std(),
                 sigma=(b.ks - n.ks.mean()) / n.ks.std())
    rows.append(r)
t = pd.DataFrame(rows)
pd.set_option("display.width", 250)
print(t.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
t.to_csv(f"{E}/chain_summary.csv", index=False)
if "V2_weighted" in perrun:
    base = perrun["V2_weighted"]
    print("\nper-run best age minus the production error model's (V2), median and range over runs:")
    for v in perrun:
        d = (perrun[v].age - base.age).dropna()
        print(f"  {v:13s} {d.median():+.2f}  ({d.min():+.1f} .. {d.max():+.1f})   best-KS change "
              f"{(perrun[v].ks - base.ks).median():+.4f}")

fig, axs = plt.subplots(1, 2, figsize=(15, 5.5), constrained_layout=True)
ax = axs[0]
y = np.arange(len(t))[::-1]
ax.barh(y, t.best_ks, color="#4a7fb5", height=0.6)
if "null_3sig" in t:
    ax.scatter(t.null_3sig, y, marker="|", s=300, color="k", zorder=3, label="its bootstrap mean + 3σ")
    ax.legend(fontsize=8, loc="lower right")
ax.set_yticks(y)
ax.set_yticklabels([f"{LABEL[v]}  ({a:.1f} Gyr, E {e:.3f})" for v, a, e in zip(t.variant, t.age, t.evi)], fontsize=8)
ax.set_xlabel("best 2D KS over 8 runs x local (age, mu, E) grid")
ax.set_title("from Martin's error model to ours, one change at a time", fontsize=10)
ax = axs[1]
for v, pr in prof.items():
    main = v.startswith("V")
    ax.plot(pr.index, pr.values, lw=2 if main else 0.8, ls="-" if main else "--", label=LABEL[v].strip())
ax.axvline(13.8, color="0.4", ls=":")
ax.set_xlabel("age (Gyr)"); ax.set_ylabel("lowest KS at that age (8 runs, any mu, E 0.050-0.100)")
ax.set_title("age profile per error model", fontsize=10); ax.legend(fontsize=7, ncol=2)
fig.savefig(f"{E}/chain_summary.png", dpi=110)
print(f"\n{E}/chain_summary.png")
