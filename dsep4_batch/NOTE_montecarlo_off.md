# Note (2026-11-18): most DSEP4 NGC 6791 test runs here did NOT read their Monte Carlo var files

DSEP4 (`~/dsep4/src`) reads the var file (fort.85, `mread`) only if `LMonteCarlo` is true, and the default is false
(`Common_variables.f90`; the call is in `parmin.f90`).

| tracks | runs | namelists | var file read |
|---|---|---|---|
| `../Tracks_dsep4` (tar) | 1 (20000) | `~/nml_dsep4` (`mcdsep/run/runmc_dsep4.sh`) | no |
| `../Tracks_dsep4_parts` | 21 | `~/nml_dsep4` (`runmc_dsep4_one.sh`) | no |
| `../Tracks_dsep4_parts_18Gyr` | 11 | most likely `~/nml_dsep4` (script not traced) | most likely no |
| `../Tracks_dsep4_mc_18Gyr` | 1 (20000, [alpha/Fe] = 0) | `~/nml_dsep4_mc` (`runmc_dsep4_mc_one.sh`) | yes |

- `~/nml_dsep4` sets `LMonteCarlo` in no namelist.
- `~/nml_dsep4_mc` sets it only in `gs98_3_DATA.nml` ([alpha/Fe] = 0).
- Runs without the read used [Fe/H], Y, [alpha/Fe] and the mixing length from the run script, and namelist defaults for
  everything else (KTTAU, diffusion, overshoot, S-factors, opacity and neutrino multipliers).

So the comparisons in `compare/` (`cmp_tracks.py`, `cmp3_mc.py`, `agree.py`), if they used these runs, mix
"Monte Carlo inputs ignored" with "DSEP4 vs dsepX".

- **Not affected:** the NGC 6791 sCMD legs. "New" is `~/dsepX_fixed` (`bd_montec.f`: `DATA LMonteCarlo/.TRUE./`),
  about 1,000 runs in `../Tracks_noco` / `../Isochrones_noco`.
- **Decision (user, 2026-11-18):** leave these tests as they are; not rerun.
- **Same bug in M15:** see `../scmd_tb/m15_dsep4_montecarlo_check`. M15 is being rerun with `LMonteCarlo` added to every
  [alpha/Fe] namelist (`../M15/dsep4/env/nml_mc`).
