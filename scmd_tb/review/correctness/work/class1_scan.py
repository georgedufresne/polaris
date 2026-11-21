import sys, os; sys.path.insert(0, ".")
from common import *
S = load("dev")
from isochronetoolbox.scmd import ast as AM
s = setup_fresh(S); ctx = S.context(s); ne = ctx.ast.nearest
runs = sorted(int(r) for r in os.listdir("/dartfs-hpc/scratch/f007058/scmd_tb/tuc_iso") if r.isdigit())
rs_pick = np.random.RandomState(0)
pick = sorted(set(rs_pick.choice(runs, 12, replace=False).tolist()) | {10099})
captured = []
orig = AM.NearestErrors.draw
def spy(self, rs, radius, v, i):
    out = orig(self, rs, radius, v, i)
    captured.append((np.array(radius), np.array(v), np.array(i), out[0]))
    return out
AM.NearestErrors.draw = spy
worst = (0, None)
for run in pick:
    for age in (8.0, 11.4, 14.0):
        try:
            _, make = ctx.model(run, age)
        except Exception as e:
            continue
        for dm, e in ((13.15, 0.05), (13.45, 0.05), (13.30, 0.0)):
            captured.clear()
            cat = make(np.random.RandomState(0), dm, e)
            if not captured:
                continue
            rad, v, i, ok = (np.concatenate([c[k] for c in captured]) for k in range(4))
            m = (v >= ne.split[0]) & (i < ne.split[1])
            if m.any():
                _, near = ne.tree.query(ne._coords(v[m], i[m], rad[m]))
                a = ctx.ast._full
                AV = a["InputF606W"].values[near] + s.zp_v
                print(f"run {run} age {age} DM {dm} E {e}: {m.sum()} class-1 stars of {len(v)}; V {v[m].min()+s.zp_v:.2f}-{v[m].max()+s.zp_v:.2f}; "
                      f"nearest-AST V median {np.median(AV):.2f}; recovered {ok[m].mean():.3f} vs {ok[~m].mean():.3f} others", flush=True)
                worst = max(worst, (m.sum(), (run, age, dm, e)))
print("max class-1 count:", worst)
