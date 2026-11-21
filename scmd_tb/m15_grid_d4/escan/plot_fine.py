#!/usr/bin/env python3
"""Fine age x E(V-I) scan of the top-20 M15 runs: where does the KS valley bottom out?"""
import glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

D = "/dartfs-hpc/scratch/f007058/scmd_tb/m15_grid_d4/escan"
NULL3 = 0.01476                                            # bootstrap mean + 3 sigma
f = pd.concat([pd.read_csv(p) for p in sorted(glob.glob(f"{D}/fine_*.csv"))], ignore_index=True)
prof = f.groupby(["run", "age", "evi"]).ks.min().reset_index()       # best mu at each (age, E)
val = prof.loc[prof.groupby(["run", "evi"]).ks.idxmin()]              # best age at each E: the valley
best = prof.loc[prof.groupby("run").ks.idxmin()].sort_values("ks")
floor = prof[prof.evi >= 0.08 - 1e-9]
bfl = floor.loc[floor.groupby("run").ks.idxmin()].set_index("run")
mu_at = f.loc[f.groupby("run").ks.idxmin()].set_index("run").mu

print("per run: free minimum vs minimum with E(V-I) >= 0.08 (fine scan, min over mu)")
print(" run    free: KS     age   E      mu   | E>=0.08: KS     age   E")
for r in best.itertuples():
    b = bfl.loc[r.run]
    print(f" {r.run}  {r.ks:.5f}  {r.age:5.1f} {r.evi:.3f} {mu_at[r.run]:.2f} | {b.ks:.5f}  {b.age:5.1f} {b.evi:.3f}")
print(f"\nall 20 runs: free minima at E {best.evi.min():.3f}-{best.evi.max():.3f}, ages {best.age.min()}-{best.age.max()};"
      f" {int((best.age >= 14.95).sum())} runs at the 15.0 Gyr edge, {int((best.evi <= 0.0301).sum())} at E 0.03")
slope = np.polyfit(val.evi, val.age, 1)[0]
print(f"valley slope d(age)/dE = {slope:.1f} Gyr/mag -> 0.01 mag of E(V-I) = {abs(slope) * 0.01:.2f} Gyr")

top = best.run.values[:6]
fig = plt.figure(figsize=(17, 9.5), constrained_layout=True)
gs = fig.add_gridspec(2, 4)
for i, run in enumerate(top):
    ax = fig.add_subplot(gs[i // 3, i % 3])
    t = prof[prof.run == run].pivot(index="age", columns="evi", values="ks")
    im = ax.pcolormesh(t.columns, t.index, t.values, cmap="Blues_r", vmin=0.015, vmax=0.08, shading="nearest")
    v = val[val.run == run].sort_values("evi")
    ax.plot(v.evi, v.age, color="0.25", lw=1, label="best age at each E")
    b = best[best.run == run].iloc[0]
    ax.plot(b.evi, b.age, "o", ms=9, mfc="none", mec="#d6604d", mew=2, label=f"minimum {b.ks:.4f}")
    ax.axvline(0.08, color="#d6604d", ls="--", lw=1.2, label="prior floor E(V-I) 0.08")
    ax.axhline(13.8, color="0.3", ls=":", lw=1.2, label="age of the Universe 13.8")
    ax.set_title(f"run {run}", fontsize=10)
    ax.set_xlabel("E(V-I)  (applied to F606W-F814W)"); ax.set_ylabel("age (Gyr)")
    if i == 0:
        ax.legend(fontsize=7, loc="lower left")
cb = fig.colorbar(im, ax=fig.axes, location="left", shrink=0.6, pad=0.01)
cb.set_label("2D KS, minimised over mu (darker = better)")

ax = fig.add_subplot(gs[0, 3])
for run, g in val.groupby("run"):
    g = g.sort_values("evi")
    ax.plot(g.evi, g.ks, color="#2166ac" if run == top[0] else "0.7", lw=2 if run == top[0] else 0.8, zorder=3 if run == top[0] else 2)
ax.axhline(NULL3, color="0.2", ls="--", lw=1); ax.text(0.031, NULL3 + 0.0006, "bootstrap null mean + 3σ", fontsize=8)
ax.axvline(0.08, color="#d6604d", ls="--", lw=1.2)
ax.set_ylim(0.0, 0.06); ax.set_xlabel("E(V-I)"); ax.set_ylabel("lowest KS over age and mu")
ax.set_title(f"valley floor vs E, 20 runs (blue: run {top[0]})", fontsize=10)

ax = fig.add_subplot(gs[1, 3])
for run, g in val.groupby("run"):
    g = g.sort_values("evi")
    ax.plot(g.evi, g.age, color="#2166ac" if run == top[0] else "0.7", lw=2 if run == top[0] else 0.8, zorder=3 if run == top[0] else 2)
ax.axvline(0.08, color="#d6604d", ls="--", lw=1.2); ax.axhline(13.8, color="0.3", ls=":", lw=1.2)
ax.set_xlabel("E(V-I)"); ax.set_ylabel("best age at that E (Gyr)")
ax.set_title(f"age-reddening degeneracy: {abs(slope) * 0.01:.2f} Gyr per 0.01 mag", fontsize=10)
fig.suptitle("M15, DSEP4: fine age x reddening scan of the 20 best runs (E(V-I) 0.030-0.110 step 0.005, ages 13.0-15.0 step 0.1, mu 15.46-15.58)")
out = f"{D}/../analysis/m15_age_reddening_scan.png"
fig.savefig(out, dpi=110); print(out)
