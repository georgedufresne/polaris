"""Check 4: CRN -- regenerate with table (rounded) values vs the grid's own floats (AST model);
crn_seed colliding with a realisation seed (Gaia fixture, cheap)."""
import sys, json, numpy as np
sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/review/statistics/dev")
from dataclasses import replace
from isochronetoolbox.scmd import Setup, context, Context, regenerate
from isochronetoolbox.scmd import bootstrap as B
from isochronetoolbox.scmd.pipeline import label
S = Setup(**{k: v for k, v in json.load(open("/dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/newmodel/setup.json")).items()
             if k in Setup.__dataclass_fields__})
print("grid floats differing from their 3-decimal label: dm", sum(float(f"{d:.3f}") != d for d in S.dm), "of", len(S.dm),
      "; ebv", sum(float(f"{e:.3f}") != e for e in S.ebv), "of", len(S.ebv))
ctx = context(S)
run, age = 10099, 11.4
_, make = ctx.model(run, age)
n_same, n_diff, ks_diff = 0, 0, []
for dm in S.dm[20:26]:
    for ebv in S.ebv[2:5]:
        g = make(np.random.RandomState(S.crn_seed), dm, ebv)                      # what grid() scores
        _, a, d, e = label(0, dm, ebv)
        r = regenerate(ctx, run, age, d, e)                                       # what the table row replays
        same = g.shape == r.shape and np.array_equal(g, r)
        n_same += same; n_diff += (not same)
        if not same:
            ks_diff.append(abs(ctx.stage3(g)["ks_statistic"] - ctx.stage3(r)["ks_statistic"]))
print(f"AST model, crn: regenerate(rounded) bit-identical to grid cell: {n_same}, different: {n_diff}; max |dKS| among different {max(ks_diff) if ks_diff else 0:.2e}")
# collision: crn_seed inside the realisation range (Gaia fixture)
FIX = "/dartfs-hpc/scratch/f007058/scmd_tb/review/statistics/dev/tests/scmd/data/iso_20001_8.3Gyr.txt"
c2 = Context(Setup(iso=FIX, dm_grid=(13.27, 13.315, 0.01), ebv_grid=(0.12, 0.165, 0.01), crn=True, crn_seed=3))
rows = B.realisations(c2, 20001, 8.3, 13.28, 0.13, range(6))
print("Gaia fixture, crn_seed=3, realisations 0-5 KS:", [round(r["ks_statistic"], 5) for r in rows])
