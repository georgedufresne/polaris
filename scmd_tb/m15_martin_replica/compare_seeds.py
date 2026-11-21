#!/usr/bin/env python3
"""Seed check for the replica: the same (run, age) units fitted twice, results/ (seed 0) and results_s1/ (seed 1:
new gp_minimize random_state, new MakeTestCMD clock seed), on runs 10000-10499. Also the 1000-run result of seed 1.

Prints: (1) per-unit scatter of his KS, DM and E between the seeds (sd of the difference / sqrt 2 = noise of one fit),
(2) best KS by age for each seed on the same runs, with each age's distance from that seed's best in units of the
per-unit noise, (3) where the overall best fit lands: seed 0, seed 1 (500 runs), seed 1 (1000 runs).
Works on partial results (paired units only). Output: logs/compare_seeds.out."""
import glob

import numpy as np
import pandas as pd

E = "/dartfs-hpc/scratch/f007058/scmd_tb/m15_martin_replica"


def load(d):
    g = pd.concat([pd.read_csv(f) for f in glob.glob(f"{E}/{d}/task_*.csv")], ignore_index=True)
    return g.drop_duplicates(["simulation_id", "age_myr"])


s0, s1 = load("results"), load("results_s1")
k = ["simulation_id", "age_myr"]
p = s0.merge(s1, on=k, suffixes=("_0", "_1"))
print(f"seed 0: {len(s0)} units; seed 1: {len(s1)} units ({s1.simulation_id.nunique()} runs); paired: {len(p)}")
near = p[np.minimum(p.ks_statistic_0, p.ks_statistic_1) < 0.025]           # units near the floor (best fits)
for name, q in (("all paired units", p), ("units with KS < 0.025", near)):
    d = q.ks_statistic_1 - q.ks_statistic_0
    print(f"  {name} ({len(q)}): KS difference median {d.median():+.5f}, sd {d.std():.5f} -> noise of one fit "
          f"{d.std() / np.sqrt(2):.5f}; DM sd/sqrt2 {(q.distance_modulus_1 - q.distance_modulus_0).std() / np.sqrt(2):.3f}, "
          f"E sd/sqrt2 {(q.extinction_vi_1 - q.extinction_vi_0).std() / np.sqrt(2):.4f}")
dn = near.ks_statistic_1 - near.ks_statistic_0
sig = dn.std() / np.sqrt(2) if len(near) > 10 else np.nan

pr = p.simulation_id.unique()
rows = {}
for tag, g in (("seed 0", s0[s0.simulation_id.isin(pr)]), ("seed 1, same runs", s1[s1.simulation_id.isin(pr)]),
               ("seed 1, 1000 runs", s1)):
    t = g.loc[g.groupby("age_gyr").ks_statistic.idxmin()].set_index("age_gyr")
    rows[tag] = t
    b = g.loc[g.ks_statistic.idxmin()]
    print(f"  {tag:18s} best KS {b.ks_statistic:.5f} at run {int(b.simulation_id)}, {b.age_gyr:.1f} Gyr, DM "
          f"{b.distance_modulus:.3f}, E {b.extinction_vi:.3f} ({g.simulation_id.nunique()} runs, {len(g)} units)")
print(f"\nbest KS by age; brackets: (KS - that column's best) / noise of one fit ({sig:.5f})")
ages = [a for a in sorted(rows["seed 0"].index) if 11.0 <= a <= 14.6]
print("  age   " + "".join(f"{t:>27s}" for t in rows))
for a in ages:
    cells = []
    for t, r in rows.items():
        if a in r.index:
            x = r.loc[a, "ks_statistic"]
            cells.append(f"{x:.5f} [{(x - r.ks_statistic.min()) / sig:4.1f}] {int(r.loc[a, 'simulation_id'])}")
        else:
            cells.append("-")
    print(f"  {a:4.1f} " + "".join(f"{c:>27s}" for c in cells))
