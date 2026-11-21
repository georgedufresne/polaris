"""Do the observed ties (magnitudes rounded to 0.001) bias any orientation?  Re-score the saved
data-vs-model pair (seed 1000) with the observed colours and F606W jittered uniformly within
their rounding interval (+-0.0005), for two jitter draws."""
import sys

import numpy as np

sys.path.insert(0, "/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox")
sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/orientation")
from isochronetoolbox.scmd.stats import rank2d  # noqa: E402
import orientlib as L  # noqa: E402

H = "/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/orientation"
z = np.load(f"{H}/data_a.jsonl.1000.npz")
ref, sim = z["ref"], z["sim"]
print("distinct colours", len(np.unique(ref[:, 0])), "distinct F606W", len(np.unique(ref[:, 1])))
base = {k: float(np.percentile(np.abs(z["d_" + k]), 95)) for k in L.KEYS}
print("as observed:", " ".join(f"{k} {v:.5f}" for k, v in base.items()))
for s in (1, 2):
    rs = np.random.RandomState(s)
    rj = ref + rs.uniform(-0.0005, 0.0005, ref.shape)
    rec, _ = L.score(rj, L.ref_ecdfs(rj, rank2d), sim, rank2d)
    print(f"jitter {s}:   ", " ".join(f"{k} {rec[k + '_p95']:.5f}" for k in L.KEYS),
          f"| p99: " + " ".join(f"{k} {rec[k + '_p99']:.5f}" for k in L.KEYS),
          f"| ying_p99 {rec['ying_p99']:.5f}")
