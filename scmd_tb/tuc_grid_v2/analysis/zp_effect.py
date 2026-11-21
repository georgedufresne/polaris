#!/usr/bin/env python3
"""Zero-point correction in the v2 grid: the map, the observed CMD before/after,
and the KS of the three best fits with and without it (grid's frozen code)."""
import json, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd
G = "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_v2"
sys.path.insert(0, f"{G}/code")
from isochronetoolbox.scmd import Setup, context, regenerate
saved = json.load(open(f"{G}/setup.json"))
s = Setup(**{k: v for k, v in saved.items() if k in Setup.__dataclass_fields__})
ctx = context(s)
ctx_raw = context(Setup(**{**{k: v for k, v in saved.items() if k in Setup.__dataclass_fields__}, "zp_correct": False}))
corr = ctx.zp_map
raw = pd.read_csv(s.obs_raw, sep=r"\s+")
x, y, v, vi = raw.x.values, raw.y.values, raw.Vvega.values, raw.VIvega.values
lo, hi = s.sgb_v - s.window_dmag, s.sgb_v + s.window_dmag
w = (v > lo) & (v < hi)
print(f"stars {len(v)}; in window {w.sum()}")
print(f"correction (window stars): median {np.median(corr[w]):+.4f}  std {np.std(corr[w]):.4f}  "
      f"16-84 {np.percentile(corr[w],16):+.4f} {np.percentile(corr[w],84):+.4f}  "
      f"1-99 {np.percentile(corr[w],1):+.4f} {np.percentile(corr[w],99):+.4f}")
print(f"fit stars: corrected {len(ctx.obs)}  uncorrected {len(ctx_raw.obs)}")
cells = [(10005, 11.1, 13.44, 0.030), (10035, 11.4, 13.38, 0.000), (10030, 12.3, 13.34, 0.005)]
rows = []
for c in cells:
    cat = regenerate(ctx, *c)
    cat_r = regenerate(ctx_raw, *c)
    rows.append(dict(run=c[0], age=c[1], mu=c[2], evi=c[3],
                     ks_corrected=ctx.stage3(cat)["ks_statistic"],
                     ks_same_scmd_vs_uncorrected_obs=ctx_raw.stage3(cat)["ks_statistic"],
                     ks_uncorrected_pipeline=ctx_raw.stage3(cat_r)["ks_statistic"]))
df = pd.DataFrame(rows); print(df.to_string(index=False))
df.to_csv(f"{G}/analysis/zp_effect.csv", index=False)
np.save(f"{G}/analysis/zp_map.npy", corr)

fig = plt.figure(figsize=(17, 5.6))
a1 = fig.add_subplot(1, 4, 1)
sc = a1.scatter(x[w], y[w], c=corr[w] * 1000, s=0.3, cmap="RdBu_r", vmin=-15, vmax=15, rasterized=True)
a1.set_aspect("equal"); a1.set_xlabel("x (px)"); a1.set_ylabel("y (px)")
a1.set_title("colour correction on the detector\n(window stars)", fontsize=10)
fig.colorbar(sc, ax=a1, shrink=0.8, label="correction (mmag), subtracted from F606W-F814W")
a2 = fig.add_subplot(1, 4, 2)
a2.hist(corr[w] * 1000, bins=np.arange(-25, 25.1, 0.5), color="#4a6fa5")
a2.set_xlabel("correction (mmag)"); a2.set_ylabel("stars")
a2.set_title(f"std {np.std(corr[w])*1000:.1f} mmag, 1-99%: {np.percentile(corr[w],1)*1000:+.0f} to {np.percentile(corr[w],99)*1000:+.0f}", fontsize=10)
for k, (lab, cc) in enumerate((("before", vi), ("after", vi - corr)), 3):
    a = fig.add_subplot(1, 4, k)
    m = w & (v > 17.4) & (v < 19.0)
    a.scatter(cc[m], v[m], s=0.2, c="k", alpha=0.3, linewidths=0, rasterized=True)
    a.set_xlim(0.50, 0.75); a.set_ylim(19.0, 17.4)
    a.set_xlabel("F606W - F814W"); a.set_title(f"upper main sequence, {lab}", fontsize=10)
    if k == 3: a.set_ylabel("F606W")
fig.tight_layout()
fig.savefig(f"{G}/analysis/zp_effect.png", dpi=120)
