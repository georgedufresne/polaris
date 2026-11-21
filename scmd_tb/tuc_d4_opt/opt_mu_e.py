#!/usr/bin/env python3
"""Continuous optimisation of the 2D KS over (distance modulus, E(V-I)) for one DSEP4 (run, age).

The grid scores the KS only on the 8 (DM, E) cells per age that stage 2's JS ranks best, on
steps of 0.01 / 0.005. M. Ying minimised the KS itself (skopt gp_minimize, 50 calls, 20 random
starts). Here, with the tuc_grid_d4 setup and its common random numbers (RandomState(crn_seed)
per evaluation, so the KS is a deterministic function of DM and E):

  grid     the grid's best cell for this (run, age)
  nm_grid  Nelder-Mead from the grid cell (the cheap candidate for running on every (run, age))
  martin   gp_minimize, 50 calls, 20 initial points, EI (his settings)
  de       differential evolution (Sobol start, population 24, up to 25 generations)
  polish   Nelder-Mead from the best point found so far
  (each method's own result and number of KS evaluations are kept, cached points included)
  seeds    the final optimum and the grid cell re-scored with seeds crn_seed+1 .. +20:
           if the optimum only beats the grid cell for crn_seed, the optimiser fitted that
           draw's noise, not the CMD

Bounds: DM 13.15-13.55, E(V-I) 0-0.06.

    opt_mu_e.py RUN AGE      -> out/opt_RUN_AGE.csv (every evaluation), out/opt_RUN_AGE.json
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution, minimize
from skopt import gp_minimize

F = os.path.dirname(os.path.abspath(__file__))
G = "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_d4"
sys.path.insert(0, os.path.join(G, "code"))
from isochronetoolbox.scmd.pipeline import Setup, context  # noqa: E402

run, age = int(sys.argv[1]), float(sys.argv[2])
s = Setup(**json.load(open(os.path.join(G, "setup.json"))))
ctx = context(s)
_, make = ctx.model(run, age)
BOUNDS = [(13.15, 13.55), (0.0, 0.06)]
evals, cache, stage = [], {}, ["grid"]


def ks(x, seed=None):
    dm, e = float(x[0]), float(x[1])
    key = (round(dm, 6), round(e, 6), seed)
    if key not in cache:
        cat = make(np.random.RandomState(s.crn_seed if seed is None else seed), dm, e)
        cache[key] = 1.0 if cat is None else float(ctx.stage3(cat)["ks_statistic"])
        evals.append(dict(stage=stage[0], dm=dm, ebv=e, seed=s.crn_seed if seed is None else seed, ks=cache[key]))
    return cache[key]


def best(stages):
    d = pd.DataFrame(evals)
    d = d[d.stage.isin(stages) & (d.seed == s.crn_seed)]
    r = d.loc[d.ks.idxmin()]
    return np.array([r.dm, r.ebv]), r.ks


t0 = time.time()
g = pd.read_csv(os.path.join(G, "combined_2d_ks_results.csv"))
g = g[(g.simulation_id == run) & np.isclose(g.age_gyr, age)]
x_grid = np.array(g.loc[g.ks_statistic.idxmin(), ["distance_modulus", "extinction_ebv"]], float)
ks(x_grid)
res = {}

stage[0] = "nm_grid"
r = minimize(ks, x_grid, method="Nelder-Mead",
             options=dict(initial_simplex=[x_grid, x_grid + [0.01, 0.0], x_grid + [0.0, 0.005]], xatol=2e-4, fatol=1e-6, maxfev=100))
res["nm_grid"] = (r.x, r.fun, r.nfev)

stage[0] = "martin"
r = gp_minimize(ks, BOUNDS, n_calls=50, n_initial_points=20, acq_func="EI", random_state=0)
res["martin"] = (np.array(r.x), r.fun, len(r.func_vals))
x_m, k_m = res["martin"][0], res["martin"][1]
t1 = time.time()

stage[0] = "de"
r = differential_evolution(ks, BOUNDS, popsize=12, maxiter=25, tol=1e-5, init="sobol", seed=1, polish=False)
res["de"] = (r.x, r.fun, r.nfev)
t2 = time.time()

stage[0] = "polish"
x0, _ = best(["grid", "nm_grid", "martin", "de"])
r = minimize(ks, x0, method="Nelder-Mead",
             options=dict(initial_simplex=[x0, x0 + [0.01, 0.0], x0 + [0.0, 0.004]], xatol=2e-4, fatol=1e-6, maxfev=100))
res["polish"] = (r.x, r.fun, r.nfev)
x_opt, k_opt = best(["grid", "nm_grid", "martin", "de", "polish"])
t3 = time.time()

stage[0] = "seeds"
sd = np.arange(1, 21) + s.crn_seed
k_opt_s = np.array([ks(x_opt, int(q)) for q in sd])
k_grid_s = np.array([ks(x_grid, int(q)) for q in sd])

pd.DataFrame(evals).assign(simulation_id=run, age_gyr=age).to_csv(os.path.join(F, "out", f"opt_{run}_{age:.1f}.csv"), index=False)
n = pd.DataFrame(evals).stage.value_counts().to_dict()
summary = dict(simulation_id=run, age_gyr=age,
               grid_dm=x_grid[0], grid_ebv=x_grid[1], grid_ks=cache[(round(x_grid[0], 6), round(x_grid[1], 6), None)],
               **{f"{m}_{q}": v for m, (x, f, nf) in res.items()
                  for q, v in (("dm", x[0]), ("ebv", x[1]), ("ks", f), ("nfev", nf))},
               opt_dm=x_opt[0], opt_ebv=x_opt[1], opt_ks=k_opt,
               seeds_opt_mean=k_opt_s.mean(), seeds_opt_std=k_opt_s.std(),
               seeds_grid_mean=k_grid_s.mean(), seeds_grid_std=k_grid_s.std(),
               n_evals=n, minutes=dict(martin=(t1 - t0) / 60, de=(t2 - t1) / 60, polish=(t3 - t2) / 60,
                                        seeds=(time.time() - t3) / 60))
json.dump(summary, open(os.path.join(F, "out", f"opt_{run}_{age:.1f}.json"), "w"), indent=1, default=float)
print(json.dumps(summary, default=float))
