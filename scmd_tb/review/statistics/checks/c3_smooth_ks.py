"""Check 3: does smoothing shift the sCMD distribution at the KS level?  Paired copy/smooth sCMDs vs a common reference."""
import sys, json, time, numpy as np, pandas as pd
sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/review/statistics/dev")
from isochronetoolbox.scmd import Setup, context
from isochronetoolbox.scmd.ast import ArtificialStars
from isochronetoolbox.scmd.pipeline import _window
from isochronetoolbox.scmd.stats import rank2d
seeds = [int(s) for s in sys.argv[1].split(",")]
refkind = sys.argv[2]                       # "copy" or "smooth"
S = Setup(**{k: v for k, v in json.load(open("/dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/newmodel/setup.json")).items()
             if k in Setup.__dataclass_fields__})
ctx = context(S)
smooth = ctx.ast
copy = ArtificialStars(smooth._full, *smooth._args, k=S.ast_k, smooth=False, metric=S.ast_metric, radial=S.radial)
run, age, dm, ebv = 10099, 11.4, 13.38, 0.015
_, make = ctx.model(run, age)
def pts(model, seed):
    ctx.ast = model
    p = ctx.points(make(np.random.RandomState(seed), dm, ebv))
    return p[_window(p, S.ks_window)]
t = time.time()
ref = pts(copy if refkind == "copy" else smooth, 3_000_000)
ref_cdf = rank2d(ref) / len(ref)
print(f"ref {len(ref)} in {time.time()-t:.0f} s", flush=True)
for s in seeds:
    kc = ctx.compare(ref, ref_cdf, np.column_stack([np.zeros(1)]) ) if False else None
    out = {}
    for name, model in (("copy", copy), ("smooth", smooth)):
        ctx.ast = model
        cat = make(np.random.RandomState(s), dm, ebv)
        out[name] = ctx.compare(ref, ref_cdf, cat)["ks_statistic"]
    print(f"RESULT ref={refkind} seed {s}: KS copy {out['copy']:.6f} smooth {out['smooth']:.6f} diff {out['smooth']-out['copy']:+.6f}", flush=True)
