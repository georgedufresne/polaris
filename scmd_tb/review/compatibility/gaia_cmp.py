import os, sys, pickle
root, out = sys.argv[1], sys.argv[2]
sys.path.insert(0, root)
import numpy as np
from isochronetoolbox.scmd import Context, Setup, process_age, regenerate
fx = os.path.join(root, "tests", "scmd", "data", "iso_20001_8.3Gyr.txt")
ctx = Context(Setup(iso=fx, dm_grid=(13.27, 13.315, 0.01), ebv_grid=(0.12, 0.165, 0.01)))
cats = {}
js, ks = process_age(ctx, 20001, 8.3, cats)
regen = {r["filename"]: regenerate(ctx, 20001, 8.3, r["distance_modulus"], r["extinction_ebv"]) for r in ks}
_, make = ctx.model(20001, 8.3)
direct = [make(np.random.RandomState(s), 13.29, 0.14) for s in (0, 5)]
pickle.dump({"js": js, "ks": ks, "cats": cats, "regen": regen, "direct": direct}, open(out, "wb"))
