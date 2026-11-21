# M15: are the F606W and F814W errors independent? (2026-11-16)

Purpose: test the advice to treat the F606W and F814W errors as independent, as M. Ying's MakeTestCMD does. The
advice says each filter is reduced on its own, so the errors are uncorrelated, and it suggests drawing each filter's
error independently from its own (radius, magnitude) distribution. Our pipeline instead takes both errors from one
artificial star (`ast_errors` "nearest"). This folder checks which model the M15 data support.

## Files

- `ast_error_correlation.py` (`logs/ast_error_correlation.out`, `ms_width.csv`) runs two checks on M. Ying's own M15
  data:
  - correlations of the two filters' artificial-star errors;
  - main-sequence colour width, observed vs predicted by each error model.
- `plot_ms_width.py` makes `ms_width.png`.

## Results

1. **Artificial stars** (`M15artstars.dat`; 83,373 of 100,000 recovered within 0.5 px and 0.75 mag): the errors
   (input - output) of the two filters are strongly correlated.
   - Overall (F606W 16.6-20.7): Spearman +0.68, Pearson +0.84 (winsorised at 5 sigma).
   - By radius: +0.85 to +0.92 in the inner third of the field (to 1,108 px), and +0.5 to +0.77 in the outer third.
     This is the signature of crowding: an artificial star has the same neighbours in both images, so a blend
     brightens it in both. Only the measurement noise, which matters most where crowding is weak, is independent.
   - Per (magnitude, radius) cell, independent draws widen the colour error by 1.3-5x. They widen it most at the
     centre: at F606W 19.5-20.0 and r < 1,108 px, the colour error is 0.0099 from one star's two errors and 0.0366 from
     independent draws.
2. **Main-sequence colour width** (raw `M15.dat`, no fiducial cut; F606W 19.4-20.6; core width 1.4826 x MAD about a
   smooth ridge). Each observed star gets its errors from artificial stars in its own (F606W, radius) cell:

| F606W | observed | one artificial star for both filters | independent draws |
|---|---|---|---|
| 19.4-19.6 | 0.0163 | 0.0089 | 0.0271 |
| 19.6-19.8 | 0.0168 | 0.0096 | 0.0302 |
| 19.8-20.0 | 0.0180 | 0.0096 | 0.0274 |
| 20.0-20.2 | 0.0188 | 0.0113 | 0.0305 |
| 20.2-20.4 | 0.0208 | 0.0119 | 0.0344 |
| 20.4-20.6 | 0.0212 | 0.0142 | 0.0402 |

   Observed width^2 = intrinsic^2 + error^2, so a model cannot predict more width than is observed.
   - Independent draws predict 1.7-1.9x the observed width at every magnitude, which excludes them for these data.
   - One star's errors for both filters leave 0.014-0.017 mag of other spread. That is about the size of the
     position-dependent offsets the zero-point map removes.
3. **Fits** (earlier runs): with Ying's isochrones, his independent error model made the best fit worse, from
   2.2 to 6.8 sigma above the null (`../m15_yingiso_err99`). It also fitted worst in `../m15_errchain`.

**Conclusion:** separate reductions make the photon and fitting noise independent between the filters. In M15, though,
crowding errors dominate, and they are shared. The artificial stars show this (correlation about +0.84), and the
observed main sequence is too narrow for independent errors. Taking both errors from one artificial star is consistent
with the data. Independent per-filter draws are not.
