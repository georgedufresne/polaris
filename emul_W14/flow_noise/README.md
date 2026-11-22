# How noisy is W v1.4, in NGC 6791 age? (2026-11-22)

**Question (the user):** "What is the noise of the emulator? It's a normalizing flow, so you can kind of do that, no?"

## Background

W v1.4 is a conditional neural spline flow (zuko NSF, 10 transforms, trained by maximum likelihood). For given physics,
EEP and initial mass, it gives a density over log Teff, log g, log L, log R, central He4 and He-core mass. Two points
about how the isochrones use it:

- **The position is the flow's mode (z = 0).** The NGC 6791 W v1.4 isochrones (`../gen_W14_iso.py`) never sample, so
  the sCMD fits never saw the flow's spread.
- **The mass comes from the age head, not the flow.** The mass at each EEP is solved from the separate age head, a
  deterministic MLP. Its quoted 90th-percentile error is 4.4 mdex in log age at the turnoff (`age_head_v6`), or about
  1 %. The flow's spread therefore covers position at fixed mass, not age-head error.

The flow was fit to deterministic DSEP tracks. Its width is the residual it could not explain (fit error plus EEP-
resampling jitter), not a Bayesian model uncertainty. The test below checks whether that width predicts the real error.

## Test (`flow_noise.py`, `flow_noise.batch`, `summarize.py`)

The test uses 36 random KTTAU 0/5 runs (`runs.txt`, seed 1122) that are in both the W14 and DSEP3 sCMD legs.

**A. The flow's spread.** At the run's 8.3 Gyr isochrone points, draw 100 realizations of each kind:
- *coherent:* one base-Gaussian draw shared by every point, so the whole curve moves together (the worst case);
- *independent:* a fresh draw at every point.

**B. The actual error.** The run's DSEP3 physics-on isochrone (`Isochrones_noco`) at 8.0, 8.3 and 8.6 Gyr.

**How each curve is fitted.** Every curve is fitted for age against the run's own central W v1.4 isochrones (7.0-9.9
Gyr in 0.05):
- the fit uses pseudo colour and magnitude, x = -5 log Teff and y = -2.5 log L;
- the box runs from 1.05 mag brighter to 0.70 mag fainter than the turnoff, like the 2D-KS box;
- points are IMF-weighted;
- distance and reddening are a free 2D shift, found by ICP;
- the age is a parabola through the rms minimum.

A self-test refits central isochrones against their own grid, and must return their own ages.

**Python.** `flow_noise.batch` uses the home conda env (`~/.conda/envs/fidanka-env`, torch 2.8.0 and zuko 1.2.0). The
scratch copy has no torch, but is otherwise identical.

## Results

See `summary.txt` (written by `summarize.py`).
