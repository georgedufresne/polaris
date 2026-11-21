"""Check 6: old settings reproduce old sCMDs bit for bit (47 Tuc preset, nearest, k=30, linear, observed, no crn)."""
import sys, numpy as np
root, out = sys.argv[1], sys.argv[2]
sys.path.insert(0, root)
from isochronetoolbox.scmd import context
from isochronetoolbox.scmd.pipeline import preset
s = preset("47tuc")
ctx = context(s)
_, make = ctx.model(10099, 11.4)
cats = {f"seed{sd}_dm{dm}_e{e}": make(np.random.RandomState(sd), dm, e) for sd, dm, e in [(7, 13.38, 0.015), (123, 13.25, 0.0)]}
it = ctx.grid(10099, 11.4)[1]
d0, e0, c0 = next(it); cats["grid_first"] = c0
np.savez(out, **cats)
print(root, {k: v.shape for k, v in cats.items()})
