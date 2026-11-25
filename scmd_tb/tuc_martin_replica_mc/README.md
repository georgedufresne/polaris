# 47 Tuc: M. Ying's 2D-KS pipeline on our Monte Carlo-on isochrones (2026-11-25)

**Purpose.** The user asked to redo our 47 Tuc analysis with Martin's pipeline and the updated isochrones. This
folder runs his own pieces end to end, as `../m15_martin_replica` did for M15. Only the isochrones are ours: the
Monte Carlo-on DSEP4 grid `../tuc_iso_d4v2_mcfix`, used by `../tuc_final_mc`. The result can be set against our
pipeline on the same isochrones (`../tuc_final_mc_refined`) and against Ying et al. 2025 (11.76 +/- 0.76 Gyr).

| piece | this replica | source |
|---|---|---|
| isochrones | ours, Monte Carlo on, in his reader's format: `iso/feh077cmd.<run>` (see below) | `../tuc_iso_d4v2_mcfix` |
| synthetic CMDs | his `MakeTestCMD.f90`, unchanged (md5 in `logs/fortran_md5.txt`, same as his file) | `Files_from_Martin/Analysis_code/cmd_fitting` |
| error model | his: each filter's completeness and error drawn separately from (radius, magnitude) bins | MakeTestCMD lines 210-240 |
| SGB alignment of the artificial-star bins | `define_errors.f90` with SGBartStar = [-13.5542, -13.2277] | his `47Tuc_AS_test.ipynb`: 17.02 - 30.5742 and 17.02 - 0.57 - 29.6777 |
| artificial-star inputs | `inputfiles/` -> `Files_from_Martin/Photometry/47Tuc_inputfiles` (Verr/Ierr bins, Distance.dat) | his |
| arguments | PDMF -0.84, binaries 0.009, 4,000,000 draws, deltaMag 3.0 | `global_var.sCMD_vars("47Tuc")`; `runsCMD.py` passes 3.0 for KS2d |
| observed stars | `47Tuc_fitstars_2dKS.dat`: 46,538 stars, 0.08 mag around his fiducial, F606W 15.02-19.02 | his |
| 2D KS | `martin_ks_verbatim.py` (copied from `../m15_martin_replica`; his `KS_2d.evaluate`, `_rank2`, `rankn` verbatim; 99th percentile of the ECDF differences) | `vorbin_chi2.py` |
| search | `gp_minimize` (EI, 50 calls, 20 initial points) over DM 13.15-13.45, E(V-I) 0-0.05 | `global_var.define_range("47Tuc")` |
| runs, ages | the 996 runs of `../tuc_final_mc`; his 41 ages, 8-16 Gyr in 0.2, below each run's age limit | |

The forced deviations are those of the M15 replica:
- `gp_minimize` gets a per-unit `random_state`;
- a non-finite KS returns 1.0;
- MakeTestCMD seeds from the clock.

## Isochrones (`write_martin_isochrones_tuc.py`, `write_iso.sbatch`, `logs/write_martin_isochrones_tuc.out`)

- **Magnitudes:** exactly our pipeline's. `read_isochrones`, `select_age(exact)` and `BCTable.absolute` from
  `../tuc_final_mc/code`, with its MIST ACS/WFC table ([Fe/H] -0.75). Written with 6 decimals; they read back
  within 5e-7 mag.
- **Output:** 995 files, 40,482 age blocks. 895 runs have all 41 ages. Run 10188 (age limit 8.0) has none, as in our
  fit.
- **Dropped points:** points outside the BC table, only the coolest upper-RGB points (e.g. run 10076 at 11 Gyr).
- **Skipped blocks** (42, `skipped_blocks.csv`):
  - 25 have fewer than 100 points (truncated isochrones);
  - 17 do not reach 3 mag above or 4 mag below his subgiant point, which his deltaMag 3.0 needs.
- **Earlier attempts:** `logs/*_try1.out` and `_try2.out` stopped on those two cases before they were handled. The
  partial files of the second attempt are in `iso_try2/`.

## Pilot (`pilot.sbatch`, job 11202786, `pilot_results/`)

| unit | KS (his statistic) | DM | E(V-I) | sCMD stars | time sCMD / rank / search |
|---|---|---|---|---|---|
| 10567, 10.2 Gyr | 0.02058 | 13.450 | 0.040 | 1,923,246 | 9 / 121 / 26 s |
| 10567, 10.4 Gyr | 0.01586 | 13.450 | 0.039 | 1,776,382 | 7 / 113 / 23 s |
| 10105, 12.2 Gyr | 0.01414 | 13.343 | 0.013 | 1,979,194 | 8 / 120 / 25 s |

- About 2.6 minutes per unit; peak memory 2.1 GB.
- Run 10567 is our pipeline's best fit (10.3 Gyr, mu 13.49). In his box it sits at the distance edge (13.45).
- His KS is not on our scale (99th percentile, 4M-draw model), so only comparisons within this run mean anything.

## Full run (`driver.sh`, log `driver.log`, PID `driver.pid`)

1. Grid: `run.sh` -> `array.sbatch` -> `replica_worker_tuc.py`, 40,482 units, from 2026-11-25 15:57.
   - Physics array 11202843 and free array 11202844.
   - Tasks that hang after "finished" (spank-auks) are cancelled.
2. `analyze_replica_tuc.py` without a null -> `logs/analyze_grid.out`.
3. Bootstrap null at the best unit (`null_worker.py`, 50 x 200 sets) -> `null/`.
4. `analyze_replica_tuc.py` -> `logs/analyze_replica_tuc.out`, `compare_pipelines.txt`.

Expected time is about 6.5 hours with all 270 workers running.
