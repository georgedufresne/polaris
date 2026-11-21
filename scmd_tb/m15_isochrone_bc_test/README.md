# M15: why Ying's isochrones differ from ours. Colour transformation or stellar models? (2026-11-17)

Purpose: Ying's published M15 isochrones and our DSEP4 rebuild share their inputs, yet they differ by the equivalent of
about +/-1 Gyr per run, and his are about 0.03 mag redder (`../m15_ying_isochrones`). Ying et al. (2025) convert with
Dartmouth PHOENIX bolometric corrections at the nearest [alpha/Fe]. Production uses one MIST ACS/WFC table at
[Fe/H] -2.25. This test converts our model points (log L, log Teff, log g) both ways and compares with his isochrones.

- PHOENIX tables: `/dartfs-hpc/rc/lab/C/ChaboyerB/Catherine/dsep3/iso/phx/Z*.HST_ACSWF` (BCv, V-F606W, V-F814W on
  Teff 2000-10000 K x log g -0.5 to 5.5; [Fe/H] -2.0 and -2.5 at [alpha/Fe] -0.2 to +0.8). In these tables
  BCv = M_V - M_bol: the Sun's entry is +0.082, giving M_V 4.83.
- `bc_test.py` (`logs/bc_test.out`, `bc_test.csv`): A = our models + MIST (production), B = our models + PHOENIX (linear
  in (Teff, log g), linear in [Fe/H] between -2.0 and -2.5 at the run's [Fe/H], nearest [alpha/Fe]), C = his isochrones.
  It uses the same 29 runs and the same feature and equivalent-age code as `../m15_ying_isochrones/compare_isochrones.py`.
  The A equivalent ages reproduce that run's exactly.

## Results

1. **The tables themselves.** PHOENIX minus MIST at [Fe/H] -2.25:
   - F606W-F814W is redder by +0.035 at the turnoff, +0.040 at M_F606W ~4.5, +0.046-0.048 at ~5.8 and the RGB base,
     and +0.053-0.056 on the lower MS.
   - F606W is fainter by +0.04-0.09.
   - Changing [alpha/Fe] from 0.2 to 0.4 changes the colour by 0.004 or less.
2. **His minus ours, median (16-84%) over 29 runs, 13.4 Gyr:**

| | A (MIST) | B (PHOENIX) |
|---|---|---|
| turnoff colour | +0.030 (+0.011..+0.049) | -0.005 (-0.025..+0.013) |
| MS colour, M_F606W 4.5 | +0.023 | -0.011 |
| MS colour, M_F606W 5.8 | +0.031 | -0.007 |
| RGB colour | +0.038 | -0.009 |
| SGB - turnoff F606W gap | -0.049 (-0.172 at 14.2 Gyr) | -0.055 (-0.177 at 14.2 Gyr) |

   The PHOENIX tables remove the systematic colour offset. The SGB-turnoff structure difference is unchanged.
3. **Our age that best matches his 13.4 Gyr isochrone:**
   - A: median 13.4, 16-84% 12.3-14.5, 17% of runs within 0.2 Gyr.
   - B: median 14.0, 12.9-14.65, 10% within 0.2 Gyr.
   - The per-run mismatch remains with matched tables, so the stellar models differ too.
4. **Calibration stars** (HIP 46120, HIP 106924; chi2 < 2 / median chi2), same 29 runs:
   - A: 59% / 0.7.
   - B: 0% / 405.
   - C (his isochrones): 21% / 210.
   - With PHOENIX, our main sequence is far too red for the parallax subdwarfs. MIST passes them.

## Conclusion

His isochrones differ from ours in two independent ways:
- **Colour transformation.** PHOENIX vs MIST explains his ~0.03 mag redder colours. The calibration stars clearly
  favour MIST (59% vs 0% pass on the same models).
- **The stellar models themselves.** The per-run turnoff and SGB structure differs even with the same tables. This
  needs his DSEP version and physics settings to diagnose.

Caveat: B uses linear interpolation on the table points. DSEP's own isochrone code may interpolate differently, which
could matter at the few-millimag level, not at the 0.03-0.05 level found here. We have not confirmed that these are the
exact tables Ying used. They are the Dartmouth PHOENIX set his paper cites.
