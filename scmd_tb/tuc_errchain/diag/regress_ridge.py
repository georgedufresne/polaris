#!/usr/bin/env python3
"""argv: code_dir setup out.npy [ridge]: the best-cell sCMD from default settings; optionally print the ridge."""
import json, sys
import numpy as np
sys.path.insert(0, sys.argv[1])
from isochronetoolbox.scmd import Setup, context
S = Setup(**{k: v for k, v in json.load(open(sys.argv[2])).items() if k in Setup.__dataclass_fields__})
C = context(S)
cat = C.model(10457, 12.6)[1](np.random.RandomState(S.crn_seed), 13.37, 0.015)
np.save(sys.argv[3], cat); print(sys.argv[1], cat.shape)
if len(sys.argv) > 4:
    r = C.rcut.ridge
    o = np.argsort(r.mag)
    print("window", C.v_window, "n_target", C.n_target)
    for k in np.linspace(0, len(o) - 1, 60).astype(int):
        print(f"  ridge colour {r.color[o][k]:.3f}  F606W {r.mag[o][k]:.3f}")
