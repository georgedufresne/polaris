# 47 Tuc: are the F606W and F814W artificial-star errors independent? (2026-11-25)

**Purpose.** The user asked whether an artificial star is really the same star in both filters, what Anderson et al.
2008 say, and what the advisor's / M. Ying's method does. This folder repeats `../m15_ast_error_correlation` for 47 Tuc
and records the literature and code checks.

## 1. What the artificial stars are (Anderson et al. 2008, AJ 135, 2055, arXiv:0804.2025)

Our files are the ACS Globular Cluster Treasury survey products described in that paper:
- `Files_from_Martin/47Tuc_data/47Tucartstars.dat` has 100,000 stars, about 10^5 per cluster as in §6.2. Its input
  magnitudes are instrumental, -5 to -17 as in §6.2.
- `47Tuc.dat` has that paper's catalogue columns (Vground, othv/othi, qfitV/qfitI, wV/wI; Table 4).
- HUGS (Piotto et al. 2015; Nardiello et al. 2018) is the later UV survey. Its catalogues also carry F606W/F814W
  from the same ACS images, but our two files have the 2008 formats above.

How each artificial star was made (§6.1-6.2, from the arXiv text):
- One input position and two input magnitudes: "a set of input parameters (x_in, y_in, mv_in, and mi_in)".
- It is added "into the raster for each exposure", in both filters, with a colour on the cluster fiducial sequence.
- Stars are inserted "one artificial star at a time", so they "can never interfere with each other".
- The patch is reduced with the same procedure as the real stars, which fits "simultaneously in all the exposures ...
  an x and a y position, and a flux in each of F606W and F814W" (§5.3).
- The routine "reports the star that was found closest to the inserted position". Recovery is "within 0.5 pixel and
  ... 0.75 magnitude".

So each row (InputX, InputY, InputF606W, InputF814W, OutputX, OutputY, OutputF606W, OutputF814W) is one star, at one
position, measured in both filters. Crowding is shared: §7.1 says a star "will be found in the same place relative to
the same neighbors in all the exposures, so any error related to that crowding will be the same for all
measurements".

The paper does not discuss the correlation between filters; the numbers below measure it.

## 2. What the advisor's / Ying's method does

**Papers.** Ying et al. 2023 (M92, arXiv:2306.02180, §III.2) and 2025 (arXiv:2505.02969, §IV) say only that errors
are "randomly selected from their observed distribution (which is a function of magnitude and distance from the
cluster center)". They do not say whether one artificial star supplies both errors.

**Code.** `Files_from_Martin/Analysis_code/cmd_fitting/MakeTestCMD.f90`, called by `runsCMD.py` / `TestCMDPAR.sh`
with 4,000,000 stars per isochrone as in the papers. For each synthetic star (lines 210-240):
1. a completeness test in the F606W (radius, F606W) bin;
2. an F606W error from that bin's histogram (`Verr???s.dat`), with its own random number;
3. a separate completeness test in the F814W (radius, F814W) bin;
4. an F814W error from that bin's histogram (`Ierr???s.dat`), with another random number;
5. a keep cut: within 0.08 mag of the isochrone in F606W (at the star's F814W) and in F814W (at its F606W).

So the advisor's code draws the two filters' errors and completeness independently. Our `ast_errors: "independent"`
reproduces this (`BandErrors.draw`, built like PhotoPropI.f90). That is the leg run in `../tuc_final_mc_indeperr`.

## 3. The data (`ast_error_correlation_tuc.py` -> `logs/ast_error_correlation_tuc.out`, `ms_width.csv`)

Zero points come from `../tuc_final_mc/setup.json` (30.5911, 29.7079). 75,275 of 100,000 stars are recovered within
0.5 px and 0.75 mag in both bands.

**Error correlation (input - output, same star).**
- Across the fit window (F606W 15.02-19.02): Spearman +0.61, winsorised Pearson +0.78.
- By radius: +0.62 to +0.83 (Pearson +0.81 to +0.88) in the inner third, +0.29 to +0.69 in the outer third.
- In the inner third below the turnoff, independent draws widen the colour error 2.8-3.3x. At F606W 18.5-19.0 the
  colour error is 0.0130 from one star's two errors and 0.0426 from independent draws.

**Blends are not the only cause.** The 30,808 best-recovered stars (output within 0.05 px of the input) still have
Spearman +0.56, Pearson +0.72. The shared position, neighbours and sky in the joint fit correlate even small errors.

| position offset (px) | n | Spearman | Pearson (w) | colour error: joint / independent |
|---|---|---|---|---|
| 0.00-0.05 | 30,808 | +0.56 | +0.72 | 0.0064 / 0.0122 |
| 0.05-0.10 | 833 | +0.62 | +0.68 | 0.031 / 0.109 |
| 0.10-0.20 | 454 | +0.68 | +0.77 | 0.041 / 0.178 |
| 0.20-0.50 | 222 | +0.77 | +0.92 | 0.042 / 0.237 |

**Main-sequence colour width.**
- Method: raw `47Tuc.dat`, no fiducial cut; core width 1.4826 x MAD about a smooth ridge. Each observed star gets
  artificial-star errors from its own (F606W 0.2 mag, radius decile) cell.

| F606W | n | observed | joint (one star) | independent | intrinsic left (joint) |
|---|---|---|---|---|---|
| 17.9-18.1 | 4,465 | 0.0141 | 0.0079 | 0.0160 | 0.0117 |
| 18.1-18.3 | 4,714 | 0.0140 | 0.0087 | 0.0194 | 0.0109 |
| 18.3-18.5 | 4,906 | 0.0149 | 0.0090 | 0.0242 | 0.0119 |
| 18.5-18.7 | 5,068 | 0.0153 | 0.0095 | 0.0256 | 0.0120 |
| 18.7-18.9 | 5,169 | 0.0166 | 0.0107 | 0.0264 | 0.0127 |
| 18.9-19.1 | 5,123 | 0.0182 | 0.0108 | 0.0288 | 0.0146 |
| 19.1-19.3 | 4,952 | 0.0189 | 0.0126 | 0.0297 | 0.0140 |
| 19.3-19.5 | 4,798 | 0.0205 | 0.0139 | 0.0311 | 0.0150 |
| 19.5-19.7 | 4,581 | 0.0215 | 0.0157 | 0.0350 | 0.0146 |
| 19.7-19.9 | 4,563 | 0.0232 | 0.0151 | 0.0378 | 0.0176 |

Observed width^2 = intrinsic^2 + error^2, so an error model cannot predict more width than is observed.
- **Independent draws predict 1.1-1.7x the observed width at every magnitude** (1.5-1.7x at F606W 18.3-19.1), so
  these data rule them out. M15 gave 1.7-1.9x.
- One star's errors for both filters leave 0.011-0.018 mag for real spread. That is plausible for 47 Tuc's multiple
  populations, small differential reddening and position-dependent zero points.
- This cell test is gentler than MakeTestCMD, whose bins are wider (10 radial x 12 magnitude bins over the whole
  field) and which bins F814W errors by F814W. Wider bins mix more radii into each error histogram, so MakeTestCMD's prediction should be wider still (not run here).

**Fits** (`../tuc_final_mc_indeperr`, 499 paired runs): the independent model's best fit is 5.9σ above its bootstrap
null, against 2.0σ for one star's errors.

## Conclusion

- An artificial star is one star, inserted at one position into both filters and measured by a joint two-filter fit
  (Anderson et al. 2008).
- Its two errors are correlated (+0.6 to +0.8; +0.56 even for the best-recovered stars).
- The advisor's MakeTestCMD draws them independently. For 47 Tuc, as for M15, that predicts a main sequence wider than
  observed and fits much worse.
- The separate per-filter reduction argument applies to photon noise. Crowding, shared neighbours and the shared
  position dominate here and are common to both filters.
