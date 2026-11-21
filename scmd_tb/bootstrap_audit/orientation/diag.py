"""Where in the CMD each orientation's large |ECDF difference| sits (saved per-point arrays)."""
import glob
import os

import numpy as np

H = "/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/orientation"
files = sorted(glob.glob(f"{H}/data_*.npz")) + sorted(glob.glob(f"{H}/grid*.npz")) + sorted(glob.glob(f"{H}/null_*.npz"))
for f in files:
    z = np.load(f)
    ref, sim = z["ref"], z["sim"]
    c, m = ref[:, 0], ref[:, 1]
    name = os.path.basename(f).replace(".jsonl", "")
    nb_r, nb_s = int((ref[:, 1] < 17.0).sum()), int((sim[:, 1] < 17.0).sum())
    print(f"\n== {name}: stars F606W < 17: reference {nb_r}, sCMD {nb_s} ({(nb_s - nb_r) / np.sqrt(nb_r):+.1f} sqrt(n_ref)); "
          f"F606W < 16: {int((ref[:, 1] < 16).sum())} vs {int((sim[:, 1] < 16).sum())}")
    print(f"   1D: colour CDF diff p95 {np.percentile(np.abs(z['d1c']), 95):.4f} max {np.abs(z['d1c']).max():.4f}; "
          f"F606W CDF diff p95 {np.percentile(np.abs(z['d1m']), 95):.4f} max {np.abs(z['d1m']).max():.4f} "
          f"(at F606W {m[np.argmax(np.abs(z['d1m']))]:.2f}, sign {np.sign(z['d1m'][np.argmax(np.abs(z['d1m']))]):+.0f})")
    for k in ("bb", "bf", "rb", "rf"):
        d = z["d_" + k]
        ad = np.abs(d)
        thr = np.percentile(ad, 95)
        top = ad >= thr
        # magnitude bands of the top-5% points
        bands = [(15.02, 16.5), (16.5, 17.3), (17.3, 18.0), (18.0, 19.02)]
        frac = " ".join(f"{np.mean((m[top] >= lo) & (m[top] < hi)):.2f}" for lo, hi in bands)
        print(f"   {k}: p95 {thr:.4f}  top-5% points: median F606W {np.median(m[top]):.2f}, median colour {np.median(c[top]):.3f}, "
              f"F606W bands [15-16.5, 16.5-17.3, 17.3-18, 18-19]: {frac}; sign(obs-model) mean {np.mean(np.sign(d[top])):+.2f}; "
              f"median F_ref there {np.median(z['F_' + k][top]):.2f}")
    # how much of each orientation's difference is the 1D marginal: rb = m1d - bb, bf = c1d - bb (exact counts)
    dx = {k: z["dx_" + k] for k in ("bb", "bf", "rb", "rf")}
    print(f"   exact-count identities (reference own ECDF has ties, so approximate): "
          f"max|dx_bb+dx_rb-d1m| {np.abs(dx['bb'] + dx['rb'] - z['d1m']).max():.4f}, "
          f"max|dx_bb+dx_bf-d1c| {np.abs(dx['bb'] + dx['bf'] - z['d1c']).max():.4f}")
