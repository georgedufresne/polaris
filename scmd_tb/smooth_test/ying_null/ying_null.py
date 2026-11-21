#!/usr/bin/env python3
"""M. Ying's resample bootstrap on our 47 Tuc data: observed stars resampled with replacement,
artificial-star errors added again (each star at its own magnitude and radius), completeness applied,
cut to the KS window, and scored against the observed stars with our 2D KS (95th percentile; exact ECDF
of the resample counted at the observed stars).  One slice of resamples per task."""
import json, sys
from dataclasses import replace
import numpy as np, pandas as pd
DEV, TASK, NTASK = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
sys.path.insert(0, DEV)
from isochronetoolbox.scmd import Setup, context
from isochronetoolbox.scmd.pipeline import _window
from isochronetoolbox.scmd.ast import centre
from isochronetoolbox.scmd.stats import rank2d
G = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02'
S0 = Setup(**{k: v for k, v in json.load(open(G + '/setup.json')).items() if k in Setup.__dataclass_fields__})
ctx = context(replace(S0, ast_errors='smooth', ast_k=100))
# the observed fit stars with positions (the same ridge selection the grid used)
raw = pd.read_csv(S0.obs_raw, sep=r'\s+')
keep = ctx.rcut.select(raw.VIvega.values, raw.Vvega.values, ctx.rcut.ridge)
d = raw[keep]
V, I = d.Vvega.values.astype(float), d.Ivega.values.astype(float)
xc, yc = centre(d.x.values.astype(np.float32), d.y.values.astype(np.float32))
R = np.hypot(d.x.values - xc, d.y.values - yc)
obs = np.column_stack([V - I, V]); w = _window(obs, S0.ks_window); obs = obs[w]; V, I, R = V[w], I[w], R[w]
assert len(obs) == len(ctx.obs3)
obs_cdf = rank2d(obs) / len(obs)
ast = ctx.ast
zpV, zpI = ast.zp['F606W'], ast.zp['F814W']
def ks_exact(model):
    allp = np.vstack([model, obs]); m = np.r_[np.ones(len(model), bool), np.zeros(len(obs), bool)]
    cdf_m = rank2d(allp, m)[~m] / len(model)
    return np.percentile(np.abs(obs_cdf - cdf_m), 95)
rows = []
configs = [('binned per-band (M. Ying)', 10), ('smoothed joint (ours)', 10), ('smoothed joint (ours)', 1)]
for k in range(TASK, 2 * NTASK, NTASK):            # two resamples per task and configuration
    for name, mult in configs:
        rs = np.random.RandomState(7_000_000 + 10 * k + mult + (0 if 'binned' in name else 3))
        idx = rs.randint(0, len(obs), mult * len(obs))
        v, i, r = V[idx], I[idx], R[idx]
        if 'binned' in name:
            okv, ev = ast.band['F606W'].draw(rs, r, v - zpV)
            oki, ei = ast.band['F814W'].draw(rs, r, i - zpI)
            ok = okv & oki
        else:
            ok, ev, ei = ast.nearest.draw(rs, r, v - zpV, i - zpI)
        v2, i2 = v - ev, i - ei                      # observed = true - (input - output), as the sCMDs
        p = np.column_stack([v2 - i2, v2])[ok]
        p = p[_window(p, S0.ks_window)]
        rows.append({'config': name, 'size': f'{mult} x N_obs', 'resample': k, 'n': len(p), 'ks': ks_exact(p)})
        print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(f'/dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/ying_null/parts/part_{TASK:03d}.csv', index=False)
