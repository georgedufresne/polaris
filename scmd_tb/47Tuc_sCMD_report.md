# 47 Tuc synthetic CMDs with artificial-star tests — what was done

*Written 2026-10-24. Code: isochrone-toolbox branch `scmd-ast`. Working files: `/dartfs-hpc/scratch/f007058/scmd_tb/`.*

## 1. What you asked for

* Add to the sCMD code what HST data need: instead of the injected empirical scatter (NGC 6791), sample each
  star's distance from the cluster centre and draw completeness and photometric errors from the
  artificial-star tests (AST). Everything after the catalogue must work exactly as for NGC 6791 (sCMD only,
  no DEBs).
* Make mock sCMDs for 47 Tuc with ~12 Gyr isochrones for runs 10000 and 10001.
* Use mu uniform 13.15–13.45 (Denissenkov et al. 2017; Brogaard et al. 2017) and E(V−I) uniform 0.00–0.05
  (Denissenkov et al. 2017).
* Run a test grid over all the tracks as for NGC 6791, and show visually what works.

## 2. What is in `Files_from_Martin`

| path | content |
|---|---|
| `47Tuc_data/47Tuc.dat` | ACS photometry, 158,156 stars: x, y, F606W (`Vvega`), F606W−F814W, F814W, errors, RA/Dec |
| `47Tuc_data/47Tucartstars.dat` | 99,999 artificial stars: input/output x, y, F606W, F814W, **instrumental** magnitudes |
| `47Tuc_data/47Tuc_data_cut.ipynb` | how the fit stars were selected (below) |
| `47Tuc_data/47Tuc_AS_test.ipynb` | AST zero points and the magnitude windows for the AST binning |
| `47Tuc_data/47Tuc_fitstars.dat` | 78,435 stars: F606W within 3.5 mag of the SGB and within 0.08 mag of the fiducial line |
| `47Tuc_data/47Tuc_fitstars_2dKS.dat` | 46,538 of them within 2 mag of the SGB (the fit sample) |
| `47Tuc_data/fiducial_lines.csv` | the fiducial line (v, vi) |
| `Photometry/*.f90` | `PhotoPropI` (bins the AST), `RadialDensity` (radial distribution) |
| `Photometry/47Tuc_inputfiles/` | their output for 47 Tuc: `VerrBoundary.dat`, `IerrBoundary.dat`, 2 × 120 sorted error files, `Distance.dat` |
| `Analysis_code/cmd_fitting/MakeTestCMD.f90` | the sCMD generator |
| `Analysis_code/cmd_fitting/vorbin_chi2.py` | M. Ying's statistics (Voronoi χ², 2D KS, KDE), `global_var.py` (per-cluster settings) |
| `M15_data`, `Photometry/M15_inputfiles` | the same for M15 (not used here) |

There is no M92 in the folder; I worked on 47 Tuc.

## 3. M. Ying's method, as read from the code

**Fit stars** (`47Tuc_data_cut.ipynb`): SGB at F606W = 17.02, F606W−F814W = 0.57 (MSTO 17.45, 0.52).
Stars within 3.5 mag of the SGB in F606W, then within 0.08 mag of the fiducial line in both bands
(|V − V_fid(I)| < 0.08 and |I − I_fid(V)| < 0.08, cubic splines of the fiducial), then within 2 mag of the
SGB → `fitstars_2dKS` (46,538).

**Artificial-star binning** (`PhotoPropI.f90`, per band): AST stars with input magnitudes inside the fit
stars' range (F606W 15.022–19.020, F814W 14.273–18.461, converted to instrumental with the zero points
30.5742 and 29.6777); a star is *recovered* if found within 0.5 px and 0.75 mag of its input (Anderson et
al. 2008); distances from the centroid of the AST stars within 200 px of (3000, 3000); 10 radial bins (the
first ends 1.5 widths out) and 12 magnitude bins, 6 brighter and 6 fainter than the short/long exposure
boundary (F606W 16.48, F814W 15.55); per bin the completeness = recovered / all and the sorted list of
(input − output) for the recovered stars.

**Radial distribution** (`RadialDensity.f90`): distances of all 78,435 fit stars from their own centroid
(stars within 200 px of 3000, 3000), sorted → `Distance.dat`.

**sCMD** (`MakeTestCMD.f90`): masses from a power-law mass function between the isochrone masses from 1 mag
fainter to 1 mag brighter than the fit region (SGB ± Δ); magnitudes by linear interpolation in mass;
optional binaries (q = 0.5–1); a distance per star from `Distance.dat` (inverse CDF); per band a completeness
test and an error from the (radius, magnitude) bin (inverse CDF of the sorted errors); then the F606W window
and a 0.08-mag cut in both bands around the *isochrone*. The AST bins are aligned to the isochrone by
matching subgiant-branch magnitudes. It writes 4 million stars per age; distance and reddening are applied
to the observations at comparison time.

**His statistics** (`vorbin_chi2.py`): Bayesian optimisation over (mu, colour excess), 2D KS at the 99th
percentile, or a KDE Jensen–Shannon distance. For 47 Tuc (`global_var.py`): [Fe/H] −0.77, binary fraction
0.009, mass-function slope −0.84, mu 13.15–13.45, colour excess 0–0.05 — the ranges you gave.

## 4. What I built (`isochronetoolbox/scmd/ast.py` + pipeline)

* `bin_artificial_stars` — `PhotoPropI.f90` in Python, in single precision with the Fortran's order of
  operations (the centroid is summed star by star, as the Fortran loop does).
* `ArtificialStars` — both bands' bins plus the radial distribution of the fit stars; `radii()` draws
  distances by inverse CDF (redrawing below the first step, as the Fortran does).
* `BandErrors.draw` — per star: radial and magnitude bin as `findIndxForErrors` assigns them, a completeness
  test, an error from the bin by inverse CDF.
* `make_scmd_ast` — one sCMD: mass function between the isochrone masses with F606W within 1 mag of the
  window, interpolation in mass, optional binaries, distance, both bands' completeness and errors, the
  F606W window (15.02 < V < 19.02, as the observed cut), the 0.08-mag cut around the sCMD's own isochrone
  with `lininterp`/`locate` ported exactly, drawn in batches until it has as many stars as the observed fit
  sample (46,538).
* `ASTContext` in `pipeline.py` — plugs this into the existing grid. Everything downstream — the (mu,
  reddening) grid in the same order, the stage-2 KDE JS ranking, the best four through the 2D KS, the
  no-survivor branch, units of work, claims, `collect`, `plot`, `validate_against_reference.py` — is the
  NGC 6791 code, unchanged.
* `--preset 47tuc` — fit window SGB ± 2 mag; mu 13.15–13.45 in 0.01 (applied to F606W) and E(V−I)
  0.000–0.050 in 0.005 (a colour shift), i.e. 31 × 11 = 341 cells per age; slope −0.84; no binaries (as in
  the NGC 6791 runs); ages 8.0–14.0 Gyr in 0.1 Gyr.
* Isochrones: MIST/pyeep files with MIST ACS/WFC bolometric corrections (tables downloaded from the new MIST
  address — fidanka's URL now returns 404; the Gaia table from the new address is byte-identical to the one
  the NGC 6791 runs used), or M. Ying's DSEP files read directly (`read_dsep`, which drops a block cut short).
* Extinction ratios, if E(B−V) is wanted instead: A_F606W = 2.88 E(B−V), A_F814W = 1.86 E(B−V) from the
  MIST tables at the 47 Tuc turnoff (Rv = 3.1; 2.84/1.85 on the RGB).

**Deliberate differences from the Fortran**

1. Magnitudes enter the AST bins through the zero points (mag_inst = mag − zp), not by matching each
   isochrone's SGB: here the grid applies mu and reddening to the sCMD, so its stars have real apparent
   magnitudes; at the right distance the two are the same.
2. Errors are applied as observed = true − (input − output). The Fortran adds (input − output), which flips
   any asymmetry of the error distribution (blends come out fainter instead of brighter). `martin_sign=True`
   restores it.
3. The Fortran `lininterp` reads `y(0)` — outside the array — for values between the first two points of a
   descending array; the port takes the first segment there. It only arises a magnitude below the fit window.

## 5. Checks

* **The binning is M. Ying's, exactly**: all 2 × 120 bins of his `47Tuc_inputfiles` — star counts,
  completeness, every error value — and the 78,435 distances of `Distance.dat` (to 5 × 10⁻⁶ px), reproduced
  from the raw AST and fit-star files (`tests/scmd/test_ast.py`). This needed his exact windows, which are the
  fit stars' magnitude range inside SGB ± 2 (from the AS-test notebook), not SGB ± 2 itself.
* `lininterp` equals a line-by-line transcription of the Fortran everywhere except the out-of-array case.
* NGC 6791 is untouched: all 7 earlier tests pass, including the 1e-12 regression pin and the KDE identity.

## 6. The first test grid: M. Ying's selection (`tuc_grid/`)

97 Monte Carlo runs (those whose tracks gave a usable isochrone set; `tuc_runs97.txt`) x 61 ages (8.0-14.0 Gyr)
x 341 (mu, E(V-I)) cells = 2.0 million sCMDs, 679 units, 139 CPU-h, ~40 min on 240 physics + 32 free CPUs,
no failures.  Observed sample: his 46,538 fit stars; sCMDs cut 0.08 mag around their own isochrone.

* Best: run 10030, 12.5 Gyr, mu 13.32, E(V-I) 0.005, KS 0.0129.  Next: 10005 (11.3 Gyr), 10080, 10091, 10033, 10070.
* Median best-fit age of the best 10-20 runs 11.4-11.7 Gyr (16-84 %: 11.0-13.1); mu 13.36-13.39; E(V-I) 0.01-0.03.
* The per-run KS valleys are broad (10.5-13.5 Gyr), and many runs' best fits sit on the edges of the priors
  (mu = 13.45 for 25 runs, E(V-I) = 0 for 32 and 0.05 for 23): the age is set mostly by the turnoff-SGB
  shape, and distance and reddening trade off against the physics.
* Fit quality against the physics is U-shaped in alpha_MLT (best ~1.4-2.0).
* The NGC 6791 3-sigma KS cut does not carry over: with N ~ 47,000 the bootstrap null is ~0.004, while the
  best fit is 0.013, so model systematics dominate; a 47 Tuc null has not been computed.

## 7. Our own fit stars: the ridge cut (`isochronetoolbox/scmd/ridge.py`)

You asked for our own fiducial through the cluster, cut around it, with the sCMDs treated the same way. The
code is your NGC 6791 ridge code, ported:

1. **Warp**: stars redder than c0 are made brighter by a smooth monotonic function of colour (PCHIP through
   knots, capped).  For 47 Tuc: your NGC 6791 knots with the colour offsets halved (47 Tuc's SGB is ~0.1
   wide in colour, NGC 6791's ~0.2), placed from each CMD's own turnoff colour (the bluest median colour of
   the 0.1-mag bins): c0 = TO + 0.04, knots at TO + 0.055/0.095/0.135/0.19 with -0.05/-0.15/-0.30/-0.45 mag,
   cap -0.5.  Without it the DP, one colour per magnitude slice, cuts the corner of the flat SGB.
2. **DP**: 120 magnitude x 240 colour cells, cost -log KDE density, max colour step 0.035 per slice, slope
   penalty 3, curvature penalty 15 (your NGC 6791 settings); no red bias.
3. **Snap**: 100 neighbours, local PCA normal, weighted median, twice.
4. **Unwarp**, then a **smoothing spline**, s = 5e-5 x points (raw CMD units), resampled at 1000 points
   uniform in arc length.
5. **Cut**: perpendicular distance, colour / 0.597 and F606W / 4.0 (the observed range in the box), kept if
   <= 0.15.  At fixed V that is +/-0.09 in colour (MS, RGB); at fixed colour +/-0.6 mag (SGB).

Where: F606W 15.02-19.02 (SGB +/- 2), colour 0.35-0.95, without the box 0.40-0.60 x 15.5-16.75 (blue
stragglers and noise blueward of the giants; as your NGC 6791 run removed the blue stragglers' corner).  The
ridge is traced over 14.72-19.32 and trimmed to the window: the DP's first and last slices sit where the
density is cut off, and the ridge hooked sideways there.

**Same treatment for the sCMDs.**  Each sCMD draws stars until it has as many in the trace box as the
observed CMD had (55,948), traces its own ridge through them with the same settings (same bandwidth, the warp
from its own turnoff), keeps the stars within 0.15 in the window and box, and is filled to the observed count
(47,413).  Before, the sCMDs were cut around their own noise-free isochrone while the observed stars were
cut around a ridge through scattered data, which is not the same selection.

**Speed.**  A ridge per sCMD means ~2 million ridges per grid.  The DP is vectorised over colour (identical to
your loop, to the bit), the snap over ridge points (identical to the scikit-learn loop to 1e-10), and the KDE
is binned (Gaussian filter of a fine histogram; its ridge agrees with scikit-learn's exact KDE to 2e-4
normalised, max 0.003; 0.06 s against 28 s).  A whole sCMD with its ridge takes ~0.25 s, as before.

**Checks.**  On NGC 6791 the port follows your saved ridge (`arc_length_continuous_ridge.csv`) along the MS,
turnoff and RGB (median offset 0.011 normalised, max 0.046), and it runs ~0.1 mag brighter along the red end
of the SGB; your saved ridge was probably traced on the full member list, which is not in your home
directory.  On 47 Tuc our ridge is within 0.012 in colour of M. Ying's fiducial everywhere and follows the
SGB stars more closely than his line.  `cut="isochrone"` still gives bit-identical sCMDs.

**What the figure shows** (`tuc_ridge_cut_final.png`): the sCMD has ~3.5x more stars outside 0.15 than the
observations (6.7 % against 1.7 %): the artificial-star error tails are heavier than the scatter of the
observed catalogue, which is probably quality-filtered while the artificial stars are not.  After the cut
both have the same number of stars, but within +/-0.09 the sCMD's distribution is still broader.

## 8. The error model and binaries

**What MakeTestCMD does** (and our `ast_errors="independent"` reproduces): per band, a completeness test and
an error drawn from that band's (radius, magnitude) bin, with separate random numbers. F606W and F814W errors
of a synthetic star therefore come from *different* artificial stars.

**Why that is wrong for colours**: one artificial star's two errors are correlated, r = 0.84 (a blend
brightens both bands). Independent draws double the colour scatter (1.48 MAD 0.0068 -> 0.0144) and put 9.6 %
of the stars more than 0.09 off in colour, against 1.7 % for real pairs. 1.7 % is what the 0.15 ridge cut
removes from the observed stars; with independent draws it removed 6.7 % from the sCMD.

**What we do now** (`ast_errors="nearest"`, the preset), as in your recipe: the 30 artificial stars nearest
in (input F606W, input F814W, radius), 0.1 mag ~ 200 px; one at random; lost if it was not recovered (0.5 px
and 0.75 mag), otherwise that star's own (dV, dI): V_obs = V_true + (V_out - V_in), I likewise, colour from
the two. For speed the 30 are those of the artificial star nearest to the synthetic star (found once).

**Binaries**: 2 % (your choice; M. Ying's global_var.py had 0.009 from Milone et al. 2012), flat q = 0.5-1,
added before the errors, as MakeTestCMD. Before this they were 0.

**What the comparisons show** (`tuc_compare_*.png`, five cells that fitted best in the first grid):

| sCMD | independent, no binaries | nearest, 2 % binaries |
|---|---|---|
| 10030, 12.5 Gyr | KS 0.0244 | KS 0.0302 |
| 10005, 11.3 Gyr | 0.0223 | 0.0286 |
| 10080, 11.2 Gyr | 0.0142 | 0.0278 |
| 10091, 11.0 Gyr | 0.0222 | 0.0238 |
| 10033, 10.8 Gyr | 0.0246 | 0.0336 |

* With independent errors the residual maps show too few sCMD stars on the ridge and too many in both wings
  (heavy tails).
* With the nearest artificial stars it flips: the sCMD core is too narrow. The observed MS is wider than
  photometric errors make it by ~0.012 in colour in quadrature (V 17.6-19: observed 0.013-0.016,
  sCMD 0.007-0.010) -- 47 Tuc's multiple populations, not in the model.
* The independent draws scored better because their extra scatter stands in for that spread. The KS
  values now carry this model-level misfit, so compare them between runs and ages, not with a null.
* Artificial stars are reused: median 3 times per sCMD, up to 18 within 500 px, where there are 2,894
  artificial stars for 8,217 cluster stars. That shows as short streaks parallel to the sequence in
  the scattered sCMD stars; it adds correlated noise, not bias.
* Runs 10030 and 10033 (12.5, 10.8 Gyr) have a giant branch too red by ~0.02; 10005 and 10091 match it.
