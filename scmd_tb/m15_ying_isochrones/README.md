# M15: Ying's published isochrones vs our DSEP4 rebuild (2026-11-11)

Purpose: compare Ying et al. (2025)'s own M15 Monte Carlo isochrones with ours for the same runs, before any
fit with them.

Data: Zenodo 10.5281/zenodo.13998637, `M15_MC_iso.tar.xz` (md5 fa0161cd4af0a189bdca232f64126011, checked).
It is unpacked to `/dartfs-hpc/scratch/f007058/M15/zenodo_ying2025/M15_MC_iso.h5` (10.2 GB).
- 10,000 sets, keyed 10000-19999 (our runs are 10000-10499).
- 41 ages, 8-16 Gyr in 0.2 Gyr steps.
- Columns: EEP, Mass, absolute F606W, and a column labelled "F814W" that holds F606W-F814W.
- The file is read with h5py (pandas "fixed" format); PyTables is not installed in fidanka-env.

## Results

1. **Same inputs.** For every run compared, his `var` table equals ours in [Fe/H], Y, [alpha/Fe], mixing
   length, both diffusion factors and envelope overshoot. The exception is core overshoot ALPHAC, which ours set
   to 0 (his values run up to 0.19).
2. **Different isochrones** (`compare_isochrones.py`, `logs/compare_isochrones.out`, `compare_isochrones.csv`;
   29 runs: 10357, 10392, 10245 and the production posterior's most-weighted runs). His minus ours, same run
   and age, median (16-84%):
   - **Turnoff colour:** +0.030 (+0.011 to +0.049) at 13.4 Gyr.
   - **Turnoff F606W:** +0.087 (+0.042 to +0.136) at 13.4 Gyr and +0.203 at 14.2 Gyr, i.e. his is fainter.
   - **Main-sequence colour:** +0.023 at F606W 4.5 and +0.031 at F606W 5.8.
   - **RGB colour:** +0.038.
3. **Equivalent age.** Which of our ages best matches his 13.4 Gyr isochrone within +/-2 mag of its SGB point,
   with free F606W and colour shifts as in the fit?
   - The median is 13.4 Gyr, but the 16-84% range is 12.3-14.5 and the full range 11.8-14.7.
   - The offset tracks the turnoff colour difference. Runs where his turnoff is 0.03-0.06 redder look like our
     11.8-12.9 Gyr; runs within 0.02 look like our 14.0-14.7.
   - No single input explains it: |Spearman| <= 0.35 for each of [Fe/H], mixing length, ALPHAC, [alpha/Fe],
     helium diffusion and Y.
   - Examples: production's best run 10357, his 13.4 = our 11.8; no-map best run 10245, his 13.4 = our 14.4.
   - So there is no overall age shift, but the ranking of runs changes. The age a fit with his isochrones would
     give cannot be predicted from this; it needs the fit.
4. **Calibration stars** (`ying_calib.py`, `logs/ying_calib.out`, `ying_calib_chi2.csv`; 500 runs, same
   recipe as `../m15_calib_stars/`):

| | chi2_min median | chi2 < 2 | p > 0.05 | colour minus HIP 46120, median (16-84%) |
|---|---|---|---|---|
| his isochrones | 76 | 19% | 24% | +0.015 (-0.001 to +0.043) |
| our isochrones | 0.7 | 59% | 67% | +0.0002 (-0.0001 to +0.013) |

   His pass rate matches what his M92 paper reports for its own sets (about 20% with chi2 < 2), which supports
   our implementation of the test. By his own test, our isochrones (with MIST bolometric corrections) match the
   parallax stars better than his do.

## What this means

His isochrones and ours share inputs but differ in shape by the equivalent of about +/-1 Gyr per run, with no
median offset. Swapping isochrones could therefore move the fitted age either way. Only the fit can tell. The
recreation run would use his isochrones, no zero-point map, his box, the 99th-percentile KS and his weighting,
with all of it set up to read his absolute magnitudes directly.

## Caveats

- The turnoff point (bluest point) is poorly defined on a flat turnoff. The equivalent-age match is the more
  robust measure.
- The match uses a colour weight of 6 (`K` in the script) and his 13.4 Gyr isochrone only.
- 29 runs for the shape comparison and 500 for the calibration test.
