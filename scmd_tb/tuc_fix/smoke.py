import sys, json, time, numpy as np
sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_fix/code")
from isochronetoolbox.scmd.pipeline import Setup, context, seed
F = "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_fix"
for v in ("A_v3", "B_zp", "C_zp_mix", "D_zp_mix_nearest"):
    s = Setup(**json.load(open(f"{F}/{v}/setup.json")))
    t = time.time(); ctx = context(s)
    age_yr, make = ctx.model(10091, 11.4)
    cat = make(np.random.RandomState(seed(10091, 11.4) if not s.crn else s.crn_seed), 13.44, 0.035)
    ks = ctx.stage3(cat)["ks_statistic"]
    print(f"{v:18s} run 10091 11.4 Gyr mu 13.44 E 0.035: KS {ks:.5f}  (n={len(cat)}, {time.time()-t:.0f} s)", flush=True)
