#!/usr/bin/env python3
"""Combine output/run_*.json from flow_noise.py into summary.txt."""
import glob
import json

import numpy as np

F = "/dartfs-hpc/scratch/f007058/emul_W14/flow_noise"
R = [json.load(open(f)) for f in sorted(glob.glob(f"{F}/output/run_*.json"))]
L = []
A = L.append
ref = R[0]["ref_age"]
A(f"W v1.4 flow spread vs actual error, NGC 6791, {len(R)} runs (KTTAU 0/5), reference age {ref} Gyr, "
  f"{R[0]['k']} draws per kind per run")
A("Age fits: turnoff box (1.05 mag brighter to 0.70 mag fainter than the turnoff), IMF-weighted, free distance + reddening "
  "shift, against the run's own central (mode) W v1.4 isochrones.")
A("")
st = np.array([[abs(v[0] - float(a)) for a, v in r["self"].items()] for r in R])
A(f"Self-test (central isochrone fitted against its own grid): max |error| {st.max():.4f} Gyr")
A("")
A("A. The flow's own spread (draws from the flow at the 8.3 Gyr isochrone's EEP/mass points)")
A(f"{'':14s}{'per-run sd (median)':>22s}{'pooled sd':>12s}{'mean - ref':>13s}{'sd log Teff @TO':>18s}{'sd log L @TO':>15s}"
  f"{'fits at grid edge':>19s}")
for kind in ("coherent", "independent"):
    ages = [np.array(r[kind]["ages"]) for r in R]
    sd = np.array([a.std(ddof=1) for a in ages]); mu = np.array([a.mean() for a in ages])
    pooled = np.concatenate([a - a.mean() for a in ages]).std()
    edge = sum(np.sum(~np.array(r[kind]["inside_grid"])) for r in R)
    sT = np.median([r[kind]["sd_logTeff_TO"] for r in R]); sL = np.median([r[kind]["sd_logL_TO"] for r in R])
    A(f"{kind:14s}{np.median(sd):>22.3f}{pooled:>12.3f}{(mu - ref).mean():>+13.3f}{sT:>18.5f}{sL:>15.4f}{edge:>19d}")
A("  coherent: one base-Gaussian draw shared by all points (whole curve moves together: the worst case);")
A("  independent: a fresh draw at every point. 'mean - ref': average fitted age of the draws minus the mode's age.")
A("")
A("B. Actual error: the run's DSEP3 physics-on isochrone (the 'New' leg) fitted with W v1.4 central isochrones")
A(f"{'DSEP3 age':>10s}{'W v1.4 age - DSEP3 age':>26s}{'sd over runs':>14s}{'SE':>8s}{'TO d log Teff (D3-W)':>23s}{'TO d log L':>12s}")
for a in ("8.0", "8.3", "8.6"):
    d = np.array([r["dsep3"][a][0] - float(a) for r in R])
    dT = np.array([r["dsep3"][a][3] for r in R]); dL = np.array([r["dsep3"][a][4] for r in R])
    A(f"{a:>10s}{d.mean():>+26.3f}{d.std(ddof=1):>14.3f}{d.std(ddof=1)/np.sqrt(len(d)):>8.3f}"
      f"{np.median(dT):>+23.5f}{np.median(dL):>+12.4f}")
A("")
d83 = np.array([r["dsep3"]["8.3"][0] - 8.3 for r in R])
coh = np.median([np.std(r["coherent"]["ages"], ddof=1) for r in R])
ind = np.median([np.std(r["independent"]["ages"], ddof=1) for r in R])
A("Comparison")
A(f"  run-to-run scatter of the actual error at 8.3 Gyr: {d83.std(ddof=1):.3f} Gyr; flow's predicted spread "
  f"{ind:.3f} (independent) to {coh:.3f} (coherent) Gyr")
A(f"  mean actual offset {d83.mean():+.3f} Gyr = {abs(d83.mean())/max(coh, 1e-9):.1f} x the coherent spread")
A("  sCMD posteriors for reference (scmd_5way_comparison.txt): W14 - New paired age difference -0.180 +/- 0.020 Gyr; "
  "per-run weighted-mean difference sd 0.35 Gyr (includes sCMD draw noise)")
open(f"{F}/summary.txt", "w").write("\n".join(L) + "\n")
print("\n".join(L))
