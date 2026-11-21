#!/usr/bin/env python3
"""Speed check for the replica's KS evaluations (one unit: run 10089, 13.4 Gyr).

His KS_2d.evaluate interpolates the sCMD ECDF (LinearNDInterpolator on ~1.7M points) at the 31,444 observed stars. scipy
locates each star by a walk through the triangulation that starts from the previous star's simplex, so the order of the
observed stars sets the cost, not the result. The 99th percentile over stars does not depend on their order either.
Check: evaluate at several (DM, E) with the observed stars in file order and sorted by (v, vi), permuting Obs_cdf the same
way; the KS values must be bit-identical. Also times skopt's own overhead (50 calls on a trivial function).
Output: logs/speed_test_<job>.out."""
import os
import subprocess
import sys
import time

import numpy as np
from skopt import gp_minimize

E = "/dartfs-hpc/scratch/f007058/scmd_tb/m15_martin_replica"
sys.path.insert(0, E)
import replica_worker as rw     # noqa: E402

wd = f"{E}/work/speed_{os.environ.get('SLURM_JOB_ID', 'local')}"
os.makedirs(wd, exist_ok=True)
if not os.path.exists(f"{wd}/inputfiles"):
    os.symlink(f"{E}/inputfiles", f"{wd}/inputfiles")
subprocess.run([f"{E}/fortran/MakeTestCMD", rw.ISO.format(10089), *rw.ARGS, "13400", rw.DMAG], cwd=wd, check=True,
               capture_output=True)
ks = rw.KS()
t = time.time()
n = ks.load(f"{wd}/mc10089.a13400")
print(f"sCMD {n} stars; rank + triangulation {time.time() - t:.0f} s", flush=True)
thetas = [(15.45, 0.085), (15.40, 0.10), (15.55, 0.12), (15.3, 0.08), (15.6, 0.15)]
t = time.time()
a = [ks.evaluate(th) for th in thetas]
ta = (time.time() - t) / len(thetas)
o = np.lexsort((ks.obs_data[:, 0], ks.obs_data[:, 1]))
ks.obs_data, ks.Obs_cdf = ks.obs_data[o], ks.Obs_cdf[o]
t = time.time()
b = [ks.evaluate(th) for th in thetas]
tb = (time.time() - t) / len(thetas)
for th, x, y in zip(thetas, a, b):
    print(f"  DM {th[0]:.2f} E {th[1]:.3f}: file order {x!r}  sorted {y!r}  identical {x == y}")
print(f"per evaluation: file order {ta:.2f} s, sorted {tb:.2f} s; all identical: {a == b}", flush=True)
t = time.time()
gp_minimize(lambda th: (th[0] - 15.45) ** 2 + (th[1] - 0.09) ** 2, [(15.3, 15.6), (0.08, 0.15)], acq_func="EI", n_calls=50,
            n_initial_points=20, verbose=False, random_state=1)
print(f"skopt overhead for 50 calls: {time.time() - t:.0f} s")
os.remove(f"{wd}/mc10089.a13400")
