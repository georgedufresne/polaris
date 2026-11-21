"""Grid cells + bootstrap rows under the OLD settings (47 Tuc preset as it ran: nearest, k 30, linear,
observed radii, no crn), with whichever tree is given; written to an npz for comparison."""
import sys, time; sys.path.insert(0, ".")
from common import *
which, cut = sys.argv[1], sys.argv[2]
S = load(which)
over = dict(ast_errors="nearest", cut=cut)
if which == "dev":
    over.update(ast_k=30, ast_metric="linear", radial="observed", crn=False)
s = setup_fresh(S, **over)
# drop the new keys for the old tree
if which == "old":
    pass
ctx = S.context(s)
out = {}
_, cells = ctx.grid(10099, 11.4)
for j, (dm, e, cat) in enumerate(cells):
    if j % 37 == 0 and cat is not None:
        out[f"grid_{j}"] = cat
    if j > 150:
        break
from isochronetoolbox.scmd import bootstrap as B
kw = dict(split_ast=False) if which == "old" else dict(resample_ast=False)
rows = B.realisations(ctx, 10099, 11.4, 13.38, 0.015, range(2), **kw)
out["boot_ks"] = np.array([r["ks_statistic"] for r in rows])
out["radii"] = ctx.ast.radii(np.random.RandomState(5), 1000) if which == "old" else ctx.ast.radii(np.random.RandomState(5), 1000, np.full(1000, 17.0))
np.savez(f"oldsettings_{which}_{cut}.npz", **out)
print(which, cut, "done", sorted(out))
