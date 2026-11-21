import argparse, os, sys, time, resource
root, d = sys.argv[1], sys.argv[2]
sys.path.insert(0, root)
from isochronetoolbox.scmd import __main__ as M
from isochronetoolbox.scmd.pipeline import context
ns = argparse.Namespace(out=d, preset=None, iso=None, obs=None, bc=None, exact_age=False, chunk=None,
                        ridge_width=None, ks_percentile=None)
s = M._setup(ns)
ctx = context(s)
t = time.time()
for i in range(3):
    r = ctx.ast.resample(1_000_000 + i)
print(f"{os.path.basename(d)}: ast_errors={s.ast_errors} k={s.ast_k} metric={s.ast_metric} radial={s.radial}: "
      f"resample {(time.time() - t) / 3:.2f} s each; n_ast {ctx.ast.n_ast}; maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.2f} GB")
