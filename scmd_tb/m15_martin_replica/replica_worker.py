#!/usr/bin/env python3
"""M15 fit with M. Ying's own pipeline, as close to his code as we can run it. One worker of a Slurm array.

Per (run, age) unit:
  1. sCMD: his MakeTestCMD.f90 (fortran/, unchanged copy, compiled with gfortran as in his TestCMDPAR.sh) on his
     isochrone (/dartfs-hpc/scratch/f007058/M15/iso_ying2025_martinfmt/feh227cmd.<run>, his Zenodo values), with his
     M15 artificial-star inputs (inputfiles/ = Files_from_Martin/Photometry/M15_inputfiles) and his M15 arguments:
     PDMF -0.99, binary fraction 0.013 (global_var.sCMD_vars), 4,000,000 draws, deltaMag 2.0 (runsCMD.py).
  2. 2D KS: his KS_2d.evaluate, _rank2 and rankn copied verbatim (martin_ks_verbatim.py) on his 2D-KS star file
     M15_fitstars_2dKS.dat (31,444 stars: 0.08 mag around his fiducial, F606W 16.64-20.64; his M15_data_cut.ipynb).
     The set-up lines of KS_2d.main are reproduced in KS.load below.
  3. Search: gp_minimize exactly as in KS_2d.main (EI, 50 calls, 20 initial points) over his M15 box
     (global_var.define_range: DM 15.3-15.6, E(V-I) 0.08-0.15).
The observed stars are reordered by (v, vi) after their ECDF is computed in file order. The KS values are
bit-identical and 77x faster (speed_test.py).
Deviations, all forced: random_state is set per unit (crc32 of "run:age") so a unit can be rerun (his is unseeded);
an evaluation that is not finite returns 1.0 instead of stopping the search; MakeTestCMD's own seed is the clock's
milliseconds, as in his code.

Output (appended unit by unit, so the array can be resubmitted): results/task_NNN.csv (best fit per unit) and
results/evals_task_NNN.csv (all 50 evaluations)."""
import argparse
import os
import subprocess
import sys
import time
import zlib

import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator
from skopt import gp_minimize

E = "/dartfs-hpc/scratch/f007058/scmd_tb/m15_martin_replica"
sys.path.insert(0, E)
import martin_ks_verbatim as mk     # noqa: E402

ISO = "/dartfs-hpc/scratch/f007058/M15/iso_ying2025_martinfmt/feh227cmd.{}"
OBS = "/dartfs-hpc/scratch/f007058/Files_from_Martin/M15_data/M15_fitstars_2dKS.dat"
ARGS = ["-0.99", "0.013", "4000000"]          # str(pdmf), str(binary), '4000000' as runsCMD.py passes them
DMAG = "2.0"
DM_MIN, DM_MAX, RED_MIN, RED_MAX = 15.3, 15.6, 0.08, 0.15
AGES = list(range(8000, 16001, 200))


class KS(mk.KS_2d):
    def __init__(self):
        self.VVI = True
        obs = pd.read_csv(OBS)[["vi", "v"]].to_numpy()          # read_obs, VVI True
        obs[:, 1] = -obs[:, 1]                                   # KS_2d.__init__: flip v
        self.obs_data = obs
        self.Obs_cdf = self.rankn(self.obs_data) / len(self.obs_data)
        # Row order only: sort the observed stars (and their ECDF values, computed above in file order) by (v, vi), so
        # scipy's point location walks between neighbouring stars. KS values are bit-identical, 77x faster
        # (speed_test.py, logs/speed_test_11066333.out).
        o = np.lexsort((self.obs_data[:, 0], self.obs_data[:, 1]))
        self.obs_data, self.Obs_cdf = self.obs_data[o], self.Obs_cdf[o]

    def load(self, path):
        """KS_2d.main up to the search."""
        sCMD = pd.read_csv(path, sep=r"\s+", names=["vi", "v"], skiprows=3)
        self.fit_values = sCMD.values
        self.fit_values[:, 1] = -self.fit_values[:, 1]
        self.ff = LinearNDInterpolator(self.fit_values, self.rankn(self.fit_values) / len(self.fit_values), fill_value=0)
        return len(self.fit_values)

    def safe(self, theta):
        v = self.evaluate(theta)
        return float(v) if np.isfinite(v) else 1.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--task", type=int, required=True)
    ap.add_argument("--ntasks", type=int, required=True)
    ap.add_argument("--units", default=None, help="run:age_myr,... (pilot); default all runs.txt x 41 ages")
    ap.add_argument("--out", default=f"{E}/results")
    ap.add_argument("--runs", default=f"{E}/runs.txt")
    ap.add_argument("--seed", type=int, default=0, help="0: random_state crc32('run:age') (first run); "
                    "else crc32('run:age:seed')")
    a = ap.parse_args()
    if a.units:
        units = [tuple(int(x) for x in u.split(":")) for u in a.units.split(",")]
    else:
        runs = [int(x) for x in open(a.runs).read().split()]
        units = [(r, g) for r in runs for g in AGES]
    mine = units[a.task::a.ntasks]
    os.makedirs(a.out, exist_ok=True)
    res_p, ev_p = f"{a.out}/task_{a.task:03d}.csv", f"{a.out}/evals_task_{a.task:03d}.csv"
    done = set()
    if os.path.exists(res_p):
        d = pd.read_csv(res_p)
        done = set(zip(d.simulation_id, d.age_myr))
    wd = f"{E}/work/{os.environ.get('SLURM_ARRAY_JOB_ID', 'local')}_{a.task:03d}"
    os.makedirs(wd, exist_ok=True)
    if not os.path.exists(f"{wd}/inputfiles"):
        os.symlink(f"{E}/inputfiles", f"{wd}/inputfiles")
    t0 = time.time()
    ks = KS()
    print(f"task {a.task}: {len(mine)} units ({len(done)} already done); observed {len(ks.obs_data)} stars; "
          f"setup {time.time() - t0:.1f} s", flush=True)
    for run, age in mine:
        if (run, age) in done:
            continue
        cage = f"{age:05d}"
        t1 = time.time()
        p = subprocess.run([f"{E}/fortran/MakeTestCMD", ISO.format(run), *ARGS, cage, DMAG], cwd=wd,
                           capture_output=True, text=True)
        out = f"{wd}/mc{run}.a{cage}"
        if p.returncode != 0 or not os.path.exists(out):
            print(f"  {run} {age}: MakeTestCMD failed ({p.returncode}): {p.stdout[-300:]} {p.stderr[-300:]}", flush=True)
            continue
        t2 = time.time()
        n = ks.load(out)
        t3 = time.time()
        res = gp_minimize(ks.safe, [(DM_MIN, DM_MAX), (RED_MIN, RED_MAX)], acq_func="EI", n_calls=50,
                          n_initial_points=20, verbose=False,
                          random_state=zlib.crc32((f"{run}:{age}" + (f":{a.seed}" if a.seed else "")).encode()))
        t4 = time.time()
        os.remove(out)
        row = pd.DataFrame([dict(simulation_id=run, age_myr=age, age_gyr=age / 1000, distance_modulus=res.x[0],
                                 extinction_vi=res.x[1], ks_statistic=res.fun, n_scmd=n, n_observed=len(ks.obs_data),
                                 t_scmd=round(t2 - t1, 1), t_rank=round(t3 - t2, 1), t_search=round(t4 - t3, 1))])
        row.to_csv(res_p, mode="a", header=not os.path.exists(res_p), index=False)
        ev = pd.DataFrame(dict(simulation_id=run, age_myr=age, distance_modulus=[x[0] for x in res.x_iters],
                               extinction_vi=[x[1] for x in res.x_iters], ks_statistic=res.func_vals))
        ev.to_csv(ev_p, mode="a", header=not os.path.exists(ev_p), index=False)
        print(f"  {run} {age}: KS {res.fun:.5f} at DM {res.x[0]:.3f}, E {res.x[1]:.3f}; sCMD {n} stars; "
              f"times {t2 - t1:.0f}/{t3 - t2:.0f}/{t4 - t3:.0f} s", flush=True)
    print(f"task {a.task} finished in {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
