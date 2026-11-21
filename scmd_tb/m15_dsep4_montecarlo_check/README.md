# Did our DSEP4 M15 runs read their Monte Carlo inputs? (2026-11-18)

Purpose: the isochrone tests showed that the whole M15 age gap with Ying et al. comes from the stellar models:
- production with Ying's isochrones gives 12.6 Gyr at null level (`../m15_final_yingiso`);
- PHOENIX tables leave our age at 14.3 (`../m15_final_phxbc`).
This checks whether our DSEP4 runs used the Monte Carlo physics in their var files.

## The finding (source and run set-up)

- **The switch:** DSEP4 (`~/dsep4/src`, the binary in `M15/dsep4/env/dsep4/bin`) reads the var file (fort.85)
  through `mread` only `if (lmontecarlo)` (`parmin.f90`). `LMonteCarlo` defaults to `.false.`
  (`Common_variables.f90:779`), and the other Monte Carlo hooks are gated the same way (`alex2`, `opal95`,
  `condopac`, `engeb`, `neut`, `rates`).
- **The namelists:** in the M15 run environment (`M15/dsep4/env/nml`, `run/subsrunmc.sh`) only
  `gs98_3_DATA.nml` ([alpha/Fe] = 0) sets `LMonteCarlo=.TRUE.`. The control namelists for [alpha/Fe] -0.2, +0.2, +0.4,
  +0.6 and +0.8 do not.
- **Consequence:** 484 of our 500 runs ([alpha/Fe] 0.2: 243, 0.4: 228, 0.6: 13) never read their var file. Only
  [Fe/H], Y, [alpha/Fe] and the mixing length, which the run script passes, vary in them. They used the physics
  namelist instead (`phys1.low.nml`): KTTAU=5, FGRY=FGRZ=1.0, the namelist S-factors, no envelope or core overshoot,
  and opacity and neutrino multipliers of 1.0.
- **Our var files** hold ALPHAC = 0 in all 500 runs. Ying's `var` tables equal ours in the other 20 inputs, but his
  ALPHAC is 0.03-0.19.

## mc_check.py (`logs/mc_check_11103532.out`, `mc_check.csv`)

1. **Var-file inputs** (median, 16-84%) vs the value the 484 runs used:
   - KTTAU: the files ask for 0 (157 runs), 1 (161) or 5 (182); all 484 used 5.
   - ALPHAE: 0.098 (0.031-0.172) vs 0.
   - FGRY, FGRZ: 0.90 (0.62-1.17) vs 1.0.
   - S-factors: pp 0.996 vs 0.968; p+C13 1.48 vs 1.27; p+O16 1.13 vs 0.99. The others are near the namelist values.
   - Opacity, triple-alpha and neutrino multipliers: about 1.0 +/- 0.2, vs 1.0.
2. **Ying's isochrones minus ours** (13.4 Gyr, ours with MIST), median (16-84%):

| | turnoff colour | SGB - turnoff F606W gap | MS colour (M_F606W 4.5) | n |
|---|---|---|---|---|
| runs that read their var file | +0.036 | -0.047 (-0.073..-0.014) | +0.024 | 16 |
| runs that did not | +0.030 | -0.050 (-0.083..-0.012) | +0.018 | 484 |

   - Among the runs that did not read, the MS colour difference depends on the KTTAU their file asked for: +0.007
     (KTTAU 0), +0.019 (1), +0.004 (5, the value they used). The surface boundary condition shifts MS colours by about
     0.015.
   - In the runs that read, the gap difference correlates only weakly with Ying's ALPHAC (Spearman +0.22).

## Conclusion

- **This is a real problem in our DSEP4 M15 grid.** 97% of runs omit the Monte Carlo physics (surface boundary
  condition, diffusion, envelope overshoot, nuclear rates, opacities). Our posterior therefore leaves most of the
  physics uncertainty out, and a third of runs used the wrong surface boundary condition. It needs fixing before the
  M15 result is final: set `LMonteCarlo=.TRUE.` in the alpha-enhanced control namelists and rerun.
- **It does not explain the difference from Ying's isochrones.** The 16 runs that read their inputs differ from his just
  as much in the turnoff/SGB structure. That points to the DSEP version and physics: our DSEP4, with Blouin20
  conductive opacities among other updates, vs his build, which his old run script suggests was dsepX-master.
