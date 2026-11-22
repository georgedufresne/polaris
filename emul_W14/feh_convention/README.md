# Does W v1.4 need [Fe/H] converted for the Z normalization? (2026-11-22)

## Background

There are two `gs98getz` builds, and Martin's gives 1.121× ours on every input:
- **Martin's** (`Catherine/run/gs98getz`) is used by the W training grid (`Catherine/run/runmc_DSEE.sh`). The training
  track headers carry Martin's Z and Y = label.
- **Ours** (`~/mcdsep/run/gs98getz`) is used by every grid we built: OG, NGC 6791 DSEP3, 47 Tuc DSEP3, and the M15 and
  47 Tuc DSEP4 grids. Their headers equal the var files exactly.

`../gen_W14_iso.py` passes the var file's [Fe/H] straight to W. So for the same var file, W modelled 12 % more Z than
DSEP3.

## Test

`feh_test.py` / `feh_test.batch` use the 36 runs of `../flow_noise/runs.txt`. For each run, two W central grids are built:
- **nominal:** [Fe/H] as written;
- **converted:** getz^-1(getz([Fe/H]) / 1.121), about -0.051 dex.

The run's DSEP3 isochrone at 8.0, 8.3 and 8.6 Gyr is fitted against each grid, using the `../flow_noise` fitter: turnoff
box, IMF weights, free distance + reddening shift.

## Result (`summary.txt`)

- **Turnoff temperature:** converting removes W's cool turnoff. DSEP3 - W goes from +0.0020 dex to -0.0006 dex, and the
  fit rms improves from 0.0023 to 0.0020 mag.
- **Isochrone-fit age:** unchanged (paired converted - nominal: -0.018 +/- 0.020 Gyr). This fit lets reddening absorb a
  uniform colour offset.
- **Why the sCMD may differ:** in the sCMD posteriors, E(B-V) differs by only 0.003 (W14 0.177, New 0.180). That is much
  less than the ~0.01 a 0.002 dex turnoff offset would need. So the sCMD fit did not simply absorb the offset in
  reddening, and the age effect has to be measured with the sCMD pipeline.

**Next (needs the user's OK):** regenerate the 36 runs' W v1.4 isochrones with converted [Fe/H], run them through the
W14 sCMD pipeline with separate output paths, and pair the result against the production W14 and DSEP3 legs.
