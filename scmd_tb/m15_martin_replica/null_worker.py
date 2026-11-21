#!/usr/bin/env python3
"""Bootstrap null for the replica, as Ying et al. (2024, Sec. 5.4; 2025, Sec. IV.2) describe it: take the isochrone with
the lowest 2D KS as the base model, simulate observed data from it on the sCMD plane, and compute his KS statistic for
each simulated data set (10,000 in all over the array).

One task: two independent sCMDs of the base unit from his MakeTestCMD (each seeds from the clock; the seeds are
checked to differ). Fake observations are drawn from sCMD A, and the expected ECDF is built from sCMD B, so the
model's own noise is in the null as it is in the fit. Each fake set holds as many stars as the observed file (31,444).
They are drawn without replacement from the sCMD-A stars inside the observed F606W range moved to the sCMD plane
(v - DM of the base fit), where the fit compares the data. Each set gets his rankn ECDF (computed before the row
reordering, as in the worker) and his evaluate at zero shift. Output: <--out, default null>/null_task_NNN.csv."""
import argparse
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd

E = "/dartfs-hpc/scratch/f007058/scmd_tb/m15_martin_replica"
sys.path.insert(0, E)
import replica_worker as rw     # noqa: E402


def scmd(wd, run, age, tag):
    cage = f"{age:05d}"
    subprocess.run([f"{E}/fortran/MakeTestCMD", rw.ISO.format(run), *rw.ARGS, cage, rw.DMAG], cwd=wd, check=True,
                   capture_output=True)
    p = f"{wd}/{tag}"
    os.replace(f"{wd}/mc{run}.a{cage}", p)
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", type=int, required=True)
    ap.add_argument("--n", type=int, required=True, help="fake data sets for this task")
    ap.add_argument("--run", type=int, required=True)
    ap.add_argument("--age", type=int, required=True, help="Myr")
    ap.add_argument("--dm", type=float, required=True)
    ap.add_argument("--out", default=f"{E}/null")
    a = ap.parse_args()
    t0 = time.time()
    wd = f"{E}/work/null_{os.environ.get('SLURM_ARRAY_JOB_ID', 'local')}_{a.task:03d}"
    os.makedirs(wd, exist_ok=True)
    if not os.path.exists(f"{wd}/inputfiles"):
        os.symlink(f"{E}/inputfiles", f"{wd}/inputfiles")
    pa = scmd(wd, a.run, a.age, "A")
    time.sleep(0.05)
    pb = scmd(wd, a.run, a.age, "B")
    A = pd.read_csv(pa, sep=r"\s+", names=["vi", "v"], skiprows=3).values
    ks = rw.KS()
    n_obs = len(ks.obs_data)
    v_obs = -ks.obs_data[:, 1]                        # observed F606W (KS stores -v)
    lo, hi = v_obs.min() - a.dm, v_obs.max() - a.dm
    nb = ks.load(pb)
    assert len(A) != nb or not np.array_equal(A[:1000], ks.fit_values[:1000] * [1, -1]), "sCMDs A and B are identical"
    pool = A[(A[:, 1] > lo) & (A[:, 1] < hi)]
    rs = np.random.RandomState(1000 + a.task)
    out = []
    for k in range(a.n):
        s = pool[rs.choice(len(pool), n_obs, replace=False)].copy()
        s[:, 1] = -s[:, 1]
        ks.obs_data = s
        ks.Obs_cdf = ks.rankn(s) / n_obs
        o = np.lexsort((s[:, 0], s[:, 1]))
        ks.obs_data, ks.Obs_cdf = s[o], ks.Obs_cdf[o]
        out.append(ks.evaluate((0.0, 0.0)))
    os.makedirs(a.out, exist_ok=True)
    pd.DataFrame(dict(task=a.task, k=np.arange(a.n), ks_statistic=out)).to_csv(f"{a.out}/null_task_{a.task:03d}.csv",
                                                                             index=False)
    for p in (pa, pb):
        os.remove(p)
    print(f"task {a.task}: base {a.run} {a.age} Myr DM {a.dm:.3f}; sCMD A {len(A)} ({len(pool)} in range), B {nb}; "
          f"{a.n} sets, KS mean {np.mean(out):.5f} sd {np.std(out):.5f}; {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
