import sys, os, numpy as np
sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/smooth_dev")
sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/smooth_dev/tests/scmd")
from isochronetoolbox.scmd import Context, Setup, regenerate, process_age
from isochronetoolbox.scmd.pipeline import label
HERE="/dartfs-hpc/scratch/f007058/scmd_tb/smooth_dev/tests/scmd"
s = Setup(iso=os.path.join(HERE,"data","iso_20001_8.3Gyr.txt"), dm_grid=(13.27,13.315,0.01), ebv_grid=(0.12,0.145,0.01), crn=True)
ctx = Context(s)
print("dm grid", repr(list(s.dm)), "ebv", repr(list(s.ebv)))
for dm, ebv, cat in ctx.grid(20001, 8.3)[1]:
    _, a, d, e = label(8.3e9, dm, ebv)
    r = regenerate(ctx, 20001, 8.3, d, e)
    k1 = ctx.stage3(cat)["ks_statistic"]; k2 = ctx.stage3(r)["ks_statistic"]
    print(repr(dm), repr(ebv), d, e, np.array_equal(cat, r), np.nanmax(np.abs(cat-r)), k1==k2)
# hand-given off-grid
r = regenerate(ctx, 20001, 8.3, 13.2804, 0.12)
c = [c for d,e,c in ctx.grid(20001,8.3)[1] if abs(d-13.28)<1e-9 and abs(e-0.12)<1e-9][0]
print("hand 13.2804: equal", np.array_equal(c, r), np.nanmax(np.abs(c-r)))
s2 = Setup(iso=s.iso, dm_grid=s.dm_grid, ebv_grid=s.ebv_grid, crn=False); ctx2=Context(s2)
r2 = regenerate(ctx2, 20001, 8.3, 13.2804, 0.12)
c2 = [c for d,e,c in ctx2.grid(20001,8.3)[1] if abs(d-13.28)<1e-9 and abs(e-0.12)<1e-9][0]
print("non-crn hand 13.2804: equal", np.array_equal(c2, r2))
