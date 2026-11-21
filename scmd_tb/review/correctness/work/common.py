import json, os, sys, time
import numpy as np, pandas as pd
R = "/dartfs-hpc/scratch/f007058/scmd_tb/review/correctness"
def load(which):
    """Import the package from dev/ or old/ (fresh interpreter per tree)."""
    sys.path.insert(0, os.path.join(R, which))
    import isochronetoolbox.scmd as S
    assert S.__file__.startswith(os.path.join(R, which)), S.__file__
    return S
def setup_fresh(S, **over):
    """The newmodel_fresh setup.json, with package-data paths re-pointed at this tree."""
    d = json.load(open("/dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/newmodel_fresh/setup.json"))
    pkg = os.path.dirname(os.path.dirname(os.path.dirname(S.__file__)))
    for k in ("completeness", "scatter_red", "scatter_blue", "bc"):
        d[k] = os.path.join(pkg, "data", *d[k].split("/data/")[1].split("/"))
    d = {k: v for k, v in d.items() if k in S.Setup.__dataclass_fields__}
    d.update(over)
    return S.Setup(**d)
