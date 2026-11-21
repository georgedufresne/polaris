# How do we match Ying et al. (2025) for M15?

You are helping a research group compare their globular-cluster age pipeline with a published one. Read
the whole brief, then answer the questions at the end. Some facts below are marked unknown. Don't fill
those in by assumption: say what each would change and how we could find out.

## Goal

Ying et al. (2025) (M. Ying, B. Chaboyer and collaborators) estimate absolute ages for eight Milky Way
globular clusters. They fit simulated colour-magnitude diagrams (sCMDs), built from Monte Carlo DSEP
isochrones, to HST photometry. We rebuilt the analysis with our own pipeline on the same photometry. For M15
(NGC 7078) we get 14.2 Gyr (68%: 14.05-14.33). The published value is about 13.2-13.4 Gyr. We want to know:

- (a) what to change to reproduce their result with the same data and the same logic;
- (b) which individual choices cause the ~0.9 Gyr difference.

M15 is the main target. 47 Tuc (NGC 104) is a second cluster where we ran the same tests. We have the
paper's sections 1-3 (introduction, data, isochrone construction) and the first part of section 4 (sCMD
construction, the Voronoi-binning and 2D KS methods, the KS bootstrap). The rest of section 4 and the
results (section 5) are not in what we have.

## Shared inputs

- **Photometry:** HST ACS Globular Cluster Treasury survey (Sarajedini et al. 2007; Anderson et al. 2008),
  F606W and F814W, VEGAmag. The M15 catalogue has 231,847 stars. Its median RA/Dec is 0.2 arcmin from the
  M15 centre.
- **Artificial-star tests (ASTs):** about 100,000 stars in instrumental magnitudes, from Ying's repository.
  The file has no RA/Dec, so its cluster identity is not independently verified. The instrumental-to-VEGA
  offsets only choose which ASTs supply a synthetic star's errors:
  - ours: 31.6176 (F606W) and 30.8912 (F814W), origin unverified;
  - Ying's: 31.638 and 30.910, matched by eye in his notebook.
- **Fiducial line and subgiant-branch magnitude:** from Ying's notebooks. V_SGB = 18.64 for M15, at the
  point on the SGB 0.05 mag redder than the turnoff.
- **Stellar-model inputs:** Ying's Monte Carlo "var" files, one set of input physics per run (Table 1 of the
  paper), and the same 39 masses from 0.2 to 3.0 Msun.

## What the paper says (sections 1-3 and part of 4)

- **Fit stars:** within +/-2 mag of V_SGB. Blue stragglers and outliers are removed by keeping stars within
  0.08 mag in F606W of the median ridgeline in an F814W-F606W magnitude-magnitude diagram.
- **Reddening:** clusters chosen for low reddening. Legnardi et al. (2023) is cited for no differential
  reddening in all eight.
- **Isochrones:** 10,000 Monte Carlo input sets per cluster. Each set gives 41 isochrones, 8-16 Gyr in
  0.2 Gyr steps, with 400 EEPs, built by the Dotter (2016) EEP method.
- **Calibration stars:** field subdwarfs with HST photometry and Gaia parallaxes (for M15: HD 46120 and
  HD 106924) are introduced "to test and calibrate the stellar isochrones". Detached eclipsing binaries are
  used for some clusters, none of them M15. For how they enter the final age, see "Final age" below.
- **Zero point:** no spatial zero-point correction is mentioned.
- **sCMD (section 4):** 4,000,000 points per isochrone, built in this order:
  - a radius from the observed radial distribution;
  - a mass from the present-day mass function;
  - binaries with the observed binary mass fraction (q flat 0.5-1, the two stars' light combined);
  - artificial-star errors and completeness.
  Distance and reddening are applied to the observed CMD to align it with the sCMD, so the sCMD is built
  once per isochrone.
- **2D KS (section 4.2):**
  - the observed ECDF is compared with the sCMD's ECDF (divide-and-conquer ranks, linearly interpolated) at
    each observed star;
  - the paper says it returns the maximum difference, but his code uses the 99th percentile;
  - (mu, E) are proposed by a Gaussian process within fixed bounds.
- **KS null (section 4.2):** the lowest-KS isochrone is the base model. Simulated observations are drawn
  from it, 10,000 resamples are drawn from the same model, and the KS is taken between the simulated
  observed and expected ECDFs. In his code the expected ECDF comes from a 1M-point sCMD and the simulated
  observations have the observed number of stars.
- **Final age:** for the Voronoi method (Ying et al. 2023), "the final age estimation will be based on the
  CMD fitting and calibration star weights". For M92, combining the fit with the calibration-star test made
  12.8 Gyr more than 10 times as likely as 14.0 Gyr. Whether the 2D KS ages are weighted the same way is
  presumably in the rest of section 4.
- **Width problem (section 4.1):** for M55 the observed CMD is wider in colour than the sCMD.
  - They rule out binaries, differential reddening (Legnardi et al. 2023) and multiple populations
    (VandenBerg et al. 2022).
  - They attribute it to PSF-model errors that artificial stars cannot capture, because the same PSF is
    used to inject and to measure them.
  - Instead of correcting it, they choose the 2D KS because it is more robust to it.
  - Our zero-point map removes the spatially coherent part of such errors. About 0.007-0.008 mag of excess
    width remains in both M15 and 47 Tuc.

## What Ying's code does (read from his repository)

**sCMD (MakeTestCMD.f90):**
- masses from a power-law mass function, magnitudes by linear interpolation in mass, binaries with
  q = 0.5-1;
- each star gets a radius from the fit stars' radial distribution (inverse CDF, unweighted);
- per band, a completeness test and an error drawn from the (radius, magnitude) AST bin, with F606W and
  F814W drawn independently;
- errors added as +(input - output), which flips any asymmetry: blends come out fainter instead of
  brighter;
- AST bins aligned to each isochrone by matching subgiant-branch magnitudes;
- an F606W window, then a 0.08 mag cut in both bands around the sCMD's own isochrone;
- 4 million stars per age.

**Settings:** M15 mass-function slope -0.99, binary fraction 0.013. 47 Tuc: -0.84 and 0.009.

**Fit (vorbin_chi2.py):** for each (run, age):
- Bayesian optimisation (gp_minimize, 50 calls) over distance modulus mu (applied to F606W) 15.3-15.6 and
  colour excess 0.08-0.15;
- the statistic is the 99th-percentile version of the 2D KS distance, in (colour, -F606W), masking zero
  CDFs, with the sCMD restricted to the observed F606W range;
- the 0.08-0.15 box is an E(V-I) range (Dutra & Bica 2000) but is applied directly as an F606W-F814W shift.
  The matching ACS colour shift is about 0.78 times that, i.e. 0.062-0.117.

**Null and uncertainty:** a 1M-star reference sCMD, and fake observations of 100k draws subsampled to the
observed number of stars.

**Observed fit stars:** his M15 files `M15_fitstars.dat` and `M15_fitstars_2dKS.dat` (31,444 stars in
F606W 16.64-20.64) are genuine M15, and his KS code reads them.

**The zero-point notebook:** his repository also has a spatial zero-point notebook (`M15_ZPCR.ipynb`) and
corrected files.
- The corrected files are M2 (NGC 7089) photometry labelled M15: their positions match an M2 catalogue
  100%, and that catalogue's median RA/Dec is 0.1 arcmin from M2.
- Neither his KS code nor the paper uses them.
- A hard-coded star count of 91,944 in an older Voronoi chi-square class matches the M2 files.

## What we do (production)

- **Fit stars:** same catalogue and window. One ridge line is traced through the data, and the same cut is
  applied to the data and to every sCMD, giving 32,085 fit stars. The group wants to keep this cut.
- **Spatial colour zero-point map:** for each star, the median offset from the ridge of its 50 nearest
  upper main-sequence stars on the detector (F606W 19.3-20.6), after subtracting the AST-predicted crowding
  shift. It is split between the bands as measured: F606W takes 0.57 of the colour shift for M15 and 0.08
  for 47 Tuc.
- **AST errors:** each synthetic star takes its (dF606W, dF814W) pair and completeness together from one of
  its 30 nearest ASTs, so the band errors stay correlated. Correct sign. Radii come from the raw photometry,
  weighted by 1/completeness. ASTs are aligned by the instrumental-to-VEGA offsets above.
- **Mass function and binaries (M15):** slope -1.00, binary fraction 0.017, from the literature.
- **Isochrones:** we rebuilt them ("DSEP4") from Ying's M15 var files, for 500 of the 10,000 runs, with ages
  in 0.1 Gyr steps and MIST ACS/WFC Vega bolometric corrections.
  - Core overshoot (ALPHAC) is set to 0. Ying's files vary it (e.g. 0.0027); this is probably irrelevant at
    M15's turnoff mass.
  - A de Boer rate multiplier line, set to 1, is added; his 21-line files don't have it.
  - Which DSEP build and which bolometric corrections produced his isochrones is unknown. We have not
    confirmed that our build reproduces his isochrones for identical inputs.
- **Fit:** the 95th-percentile 2D KS. A full (mu, E) grid (KDE pre-ranking, then KS), then a KS refinement
  by 8-neighbour descent with E 0.000-0.120 in 0.001 steps.
- **Uncertainty:** a 10,000-realisation bootstrap with fresh reference sCMDs at the best cell. The posterior
  keeps (run, age) cells with KS at or below the null mean + 3 sigma, weighted by the null density at their
  KS. No calibration-star or eclipsing-binary weighting.

**Results (M15):**
- best KS 0.0137, 2.3 sigma above the bootstrap null (0.0076 +/- 0.0026);
- age 14.18 (14.03-14.32) Gyr with no reddening box;
- 14.20 with mu 15.3-15.6 and E 0.062-0.117 (the converted box);
- 14.24 with Ying's box unconverted (E 0.08-0.15);
- [Fe/H] -2.31, E(F606W-F814W) 0.076, mu 15.515.

**47 Tuc:** 11.37 Gyr with priors.

**Known degeneracy:** within a run, 0.01 mag of E trades for about 0.27 Gyr of age; lower E gives an older
age.

## Tests done (paired, one change at a time, same runs)

| change | effect on M15 age | fit quality |
|---|---|---|
| Ying's AST error model (per-band independent errors, completeness product, unweighted radii) vs ours, 8 runs | per-run best age +0.35 Gyr (older); 47 Tuc -0.5 Gyr | worse: 6.4 vs 3.7 sigma above the null (47 Tuc 6.7 vs 5.2) |
| zero-point map on vs off, 8 runs | +0.05 Gyr; 47 Tuc -0.25 Gyr | better: 3.5 -> 3.0 sigma (47 Tuc 5.0 -> 3.0) |
| Ying's zero-point recipe vs ours, all 500 runs, full pipeline | +0.03-0.04 Gyr | best KS 0.0129 vs 0.0137 |
| Ying's reddening box (unconverted) vs none | 14.24 vs 14.18 | - |
| Ying's box and no zero-point map, all 500 runs, full pipeline (2026-11-11) | best fit still 14.0 Gyr, but the KS-age curve flattens: 12.3 Gyr only 0.0012 worse than 14.0 (0.008 with the map). No cell inside 3 sigma; lenient posteriors 14.11 (13.08-14.34) and 13.11 (12.11-14.12) | worse: best KS 0.0190, 4.3 sigma (0.0137, 2.3 sigma with the map) |

The 8-run map test above used a local grid (12.6-14.8 Gyr, E 0.05-0.10). The full-grid run supersedes it:
without the map, 14.0-14.5 Gyr fits get worse by 0.003-0.006 in KS while 11-13 Gyr fits stay the same or
improve slightly. The young fits sit at mu 15.60 (top of Ying's box) and E 0.080 (bottom).

**The spatial map is real.** Its offsets have rms 0.008 mag: 4 times the noise in M15 and 7 times in
47 Tuc.

**It is not reddening.** In each patch of the detector, the subgiant-branch magnitude offset was compared
with the main-sequence colour offset:
- M15: slope 5.34 +/- 0.32, against 9.46 predicted for reddening;
- 47 Tuc: slope 2.54 +/- 0.19, against 6.94.

## Not yet tested

- The 2D KS percentile: 99th (Ying) vs 95th (ours).
- Ying's two different cuts (data around the data's ridgeline, sCMD around its own isochrone) vs our single
  shared cut.
- The isochrone build: his original isochrones vs our DSEP4 rebuild for the same var files; 500 runs vs
  10,000.
- The search: 50-call Bayesian optimisation inside the box vs a full grid plus refinement.
- His error sign and SGB-matched AST alignment. Whether our 8-run test of his error model included these is
  unconfirmed.
- sCMD size. Every sCMD of ours has exactly as many stars as the observed sample (32,085 for M15), and our
  null compares N stars with N stars. His fit sCMD has 4M points, and his null compares N stars with 1M.
  Our model ECDF therefore carries sampling noise about as large as the data's. Our null KS std (0.0026)
  is larger than the 0.0012 that separates 14.0 from 12.3 Gyr without the zero-point map.
- Mass function and binaries (small differences for M15).
- How the paper turns per-run fits into an age. Calibration-star weights are confirmed for the Voronoi
  method and likely for the 2D KS. Unknown: how runs are aggregated, and whether there are distance or
  reddening priors (the rest of section 4 and section 5, which we don't have). Our posterior has no
  calibration-star weighting. Table 2 of the paper gives the M15 calibration stars' F606W and F606W-F814W:
  HD 46120 5.787(3), 0.566(2); HD 106924 6.041(4), 0.601(5).

## Constraints

- **Compute:** Slurm, up to 272 CPUs at once. A full 500-run grid plus refinement takes about 9 hours; an
  8-run local test takes well under an hour.
- **Test design:** we prefer paired, one-change tests against production on the same runs, so that each
  difference is measured on its own. A separate full-replication track (Ying's pipeline end to end) is
  acceptable if you think it is needed.

## Questions

1. Rank the untested differences by how likely each is to move the M15 age by ~0.9 Gyr toward 13.3 Gyr.
   Give the expected sign, rough size and your reasoning. The error-model and zero-point changes we tested
   either do nothing or push the age older.
2. Propose an ordered test plan to reproduce Ying et al.'s M15 age with the same data and logic. For each
   step, say:
   - what to change;
   - how to validate it (for example, reproduce one of his per-(run, age) best fits);
   - whether an 8-run local test is enough or the full 500 runs are needed.
3. Should we first reproduce his result end to end, using his code paths and keeping choices that look
   like errors, and only then switch his choices to ours one at a time? Or the reverse? Explain the
   trade-off.
4. Which of his choices look like errors rather than defensible choices: the error sign, the unconverted
   reddening box, different cuts for data and model, the 91,944 constant, the zero-point notebook run on M2
   data, the SGB-matched AST alignment? For each, should we replicate it to match his number and report the
   corrected version separately?
5. What exactly should we ask M. Ying (files, settings, the section 4-5 procedure) to settle the remaining
   unknowns?
6. If the paper weights Monte Carlo runs by calibration stars or uses a distance prior, how could each shift
   the M15 age relative to our unweighted posterior? How would you implement the same weighting in our
   pipeline?

Answer in that order. Keep claims about the paper to what is stated above, and mark anything else as an
assumption.
