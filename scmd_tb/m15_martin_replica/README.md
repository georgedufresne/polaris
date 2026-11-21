# M15: Ying et al. (2025)'s 2D-KS fit run with his own code (2026-11-14)

Purpose: rerun M. Ying's M15 2D-KS pipeline as closely as his code allows. With his isochrones, our pipeline gives
12.8 Gyr (`../m15_yingiso_err99`), about 0.6 Gyr below his published 13.39 +/- 0.39. This folder runs his own pieces
end to end. If it reproduces 13.4, his pieces are swapped back for ours one at a time to find which one moves the age.
If it does not, the gap lies outside the code we have.

| piece | this replica | source |
|---|---|---|
| isochrones | his, `M15/iso_ying2025_martinfmt/feh227cmd.<run>` (written by `write_martin_isochrones.py` in his reader's fixed format) | Zenodo 10.5281/zenodo.13998637 |
| synthetic CMDs | his `MakeTestCMD.f90`, unchanged (`fortran/`, md5 in `logs/fortran_md5.txt`), compiled as in his `TestCMDPAR.sh`: 4,000,000 draws, PDMF -0.99, binaries 0.013, deltaMag 2.0. It includes his subgiant-branch alignment of the artificial-star bins and his 0.08 mag cut around each isochrone. | `cmd_fitting/`, `global_var.sCMD_vars("M15")`, `runsCMD.py` |
| artificial-star inputs | `inputfiles/` = `Files_from_Martin/Photometry/M15_inputfiles` (Verr/Ierr bins, Distance.dat) | his |
| observed stars | `M15_fitstars_2dKS.dat`: 31,444 stars, 0.08 mag around his fiducial, F606W 16.64-20.64 | his `M15_data_cut.ipynb` |
| 2D KS | `KS_2d.evaluate`, `_rank2`, `rankn` copied verbatim into `martin_ks_verbatim.py` (lines 30-88 and 589-606 of `vorbin_chi2.py`); the set-up lines of `KS_2d.main` are in `replica_worker.KS.load` | `vorbin_chi2.py` |
| search | `gp_minimize` as in `KS_2d.main` (EI, 50 calls, 20 initial points) over DM 15.3-15.6, E(V-I) 0.08-0.15 | `global_var.define_range("M15")` |
| runs, ages | the same 500 runs as our M15 runs; his 41 ages, 8-16 Gyr in 0.2 | |

Forced deviations: `gp_minimize` gets a per-unit `random_state` (crc32 of "run:age") so a unit can be rerun (his is
unseeded). A non-finite KS value returns 1.0 instead of stopping the search. MakeTestCMD seeds from the clock, as in
his code.

## Files

- `write_martin_isochrones.py` (`logs/write_martin_isochrones.out`): 500 files, 20,500 age blocks. Every header parses
  back with MakeTestCMD's columns, and every row matches the source.
- `martin_ks_verbatim.py`: generated from `vorbin_chi2.py` (the copied line ranges are asserted when it is written).
- `replica_worker.py`: one Slurm task (MakeTestCMD -> his KS -> gp_minimize per unit). It writes
  `results/task_NNN.csv` (best fit per unit) and `results/evals_task_NNN.csv` (all 50 evaluations), unit by unit.
- `pilot.sbatch` (`logs/pilot_11065857.out`, `pilot_results/`), `speed_test.py`, `array.sbatch`, `run.sh`.

## Pilot (job 11065857, three units)

| unit | KS (his statistic) | DM | E(V-I) | sCMD stars | time sCMD / rank / search |
|---|---|---|---|---|---|
| 10089, 12.6 Gyr | 0.01549 | 15.480 | 0.090 | 1,697,808 | 8 / 104 / 719 s |
| 10089, 13.4 Gyr | 0.01816 | 15.406 | 0.080 | 1,691,793 | 8 / 103 / 654 s |
| 10027, 12.8 Gyr | 0.01988 | 15.406 | 0.080 | 1,671,376 | 7 / 102 / 600 s |

Peak memory 1.9 GB. His KS is not on the same scale as ours (one ECDF orientation, a 4M-draw model), so only
comparisons within this run are meaningful.

## Speed check (`speed_test.py`, `logs/speed_test_11066333.out`)

Sorting the observed stars by (v, vi) leaves the KS bit-identical at five (DM, E) points, from (15.30, 0.080) to
(15.60, 0.150). One evaluation takes 0.14 s instead of 10.8 s, because scipy starts each point's search where the last
one ended. skopt's own overhead is 5 s per 50 calls. One unit takes about 2 minutes (sCMD 8 s, rank and triangulation
95 s, search about 12 s). The worker computes the observed ECDF in file order, then reorders.

## Full run

Submitted 2026-11-14 16:42 (`run.sh`, `run.log`): physics array 11066361 (240 tasks) and free array 11066362 (32 tasks),
20,500 units. Expected time is about 2.5 hours.

## Results (2026-11-15)

- **Grid:** finished 2026-11-14 at about 18:52 (20,500 units, no duplicates, no errors). 86 tasks stayed RUNNING after
  writing "finished": Slurm's spank-auks plugin crashed at exit (`munmap_chunk(): invalid pointer`), and no user
  processes were left on the nodes. Every one had written all its units, so they were cancelled on 2026-11-15
  (`run.log`).
- **Null** (`null_worker.py`, array 11079489; base unit 10411, 13.2 Gyr, DM 15.429): 10,000 sets, KS
  0.00442 +/- 0.00115, maximum 0.0107, 3 sigma line 0.0079.
- **Best fit:** KS 0.01523 at run 10411, 13.2 Gyr, DM 15.429, E 0.080 (the floor of his box), 9.4 sigma above the
  null. No run is within 3 sigma.
- **Best KS by age** (`logs/analyze_replica.out`; brackets: (KS - best) / null sd):

| age (Gyr) | 11.2 | 11.4 | 11.8 | 12.0 | 12.4 | 12.8 | 13.0 | 13.2 | 13.4 | 13.6 | 14.0 | 14.4 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| best KS | 0.0168 | 0.0171 | 0.0164 | 0.0159 | 0.0157 | 0.0155 | 0.0155 | 0.0152 | 0.0157 | 0.0164 | 0.0168 | 0.0182 |
| (KS - best) / null sd | 1.4 | 1.7 | 1.0 | 0.5 | 0.4 | 0.2 | 0.2 | 0 | 0.4 | 1.0 | 1.4 | 2.6 |

- **Posterior** (rules of Ying 2023/2024/2025; the papers give no formula):
  - No cut, per isochrone or per run: KDE density gives 13.20 (all weight on run 10411). The fitted-normal density gives
    13.11-13.14 +/- 0.2 (1.4-1.7 effective runs). Calibration weights change neither.
  - Survival fraction and the M92 3-sigma rule: zero weight everywhere.
- **Comparison:**

| | age (Gyr) | DM | E(V-I) |
|---|---|---|---|
| Ying et al. 2025 Table 7, 2D KS | 13.39 +/- 0.39 | 15.49 +/- 0.03 | 0.08 +/- 0.01 |
| this replica (his code and inputs, 500 of his runs) | 13.2 | 15.43 | 0.080 |
| `../m15_yingiso_err99` (our pipeline, his isochrones, his errors, 99th-percentile KS) | 12.8 | 15.47 | 0.081 |

- **Conclusion:**
  - His own code on his inputs reproduces his published age, within half his quoted uncertainty, and his E.
  - Our pipeline with his isochrones gives 12.8. In both pipelines the best fits from 12.0 to 13.4 Gyr lie within
    0.5 null sd of each other, so the 0.4 Gyr between 12.8 and 13.2 is at the level of picking the minimum in a flat
    valley. His pipeline penalises 11.x a little more and 13.4 a little less.
  - In his pipeline too, every fit is far outside his null. His weighting then concentrates on one or two
    isochrones, so the code here cannot reproduce his +/- 0.39 Gyr. His null construction or his other 9,500 runs
    would have to supply it.

## Repeat: 1000 runs, different seed (submitted 2026-11-15 12:11)

The user asked to repeat with 1000 runs and a different seed. `runs_1000.txt` holds the first 500 runs again plus
10500-10999. That gives a seed check (same 500 runs) and a run-count check (500 more) in one grid.

- `write_more_isochrones.py` (`logs/write_more_isochrones.out`) writes runs 10500-10999 straight from his HDF5 with
  the same writer. Its output for the first 500 is byte-identical to the files already used. It also writes
  `ying_calib_chi2_1000.csv` (calibration weights for all 1000 runs), which matches `../m15_ying_isochrones` to
  2e-13 over the first 500.
- `replica_worker.py` has `--runs` and `--seed`. Seed 0 (the default) keeps the first run's
  random_state crc32("run:age"). Seed 1 uses crc32("run:age:1"). MakeTestCMD reseeds from the clock as always.
- `run_s1.sh` (`run.log`): physics array 11079549, free array 11079550, 41,000 units, output in `results_s1/`.
- `analyze_replica.py` takes `--results/--null/--cal/--tag`; its defaults reproduce `logs/analyze_replica.out`.
  `compare_seeds.py` (`logs/compare_seeds.out`) pairs the two seeds unit by unit.

### Results of the repeat (2026-11-16)

- **Grid:** 41,000 units, no errors. 79 tasks hung after finishing (spank-auks crash at exit). Every one had written
  all its units, so they were cancelled.
- **Null** (`null_s1/`, array 11082171; base unit 10336, 13.2 Gyr, DM 15.443): 0.00442 +/- 0.00111, the same as
  the first run's (0.00442 +/- 0.00115).
- **Seed check** (`logs/compare_seeds.out`): 20,500 paired units.
  - For units near the floor (KS < 0.025), refitting changes the KS by 0.00107 per fit, about the null's spread.
    DM changes by 0.007 and E by 0.0004.
  - The best fit moves to a different isochrone (seed 0: run 10411; seed 1: run 10336) but stays at 13.2 Gyr. The 500
    new runs (10500-10999) do not beat it.
- **Best KS by age** ((KS - best) / noise of one fit):

| age (Gyr) | 11.2 | 11.6 | 12.0 | 12.4 | 12.8 | 13.0 | 13.2 | 13.4 | 13.6 | 14.0 |
|---|---|---|---|---|---|---|---|---|---|---|
| seed 0, 500 runs | 1.5 | 1.7 | 0.6 | 0.4 | 0.2 | 0.2 | 0 | 0.4 | 1.1 | 1.5 |
| seed 1, 1000 runs | 1.5 | 1.3 | 0.7 | 0.5 | 0.5 | 0.3 | 0 | 0.2 | 0.9 | 1.8 |

- **Posterior** (`logs/analyze_replica_s1.out`, 1000 runs):
  - KDE reading: 13.20 Gyr (all weight on run 10336, DM 15.443, E 0.085).
  - Fitted-normal reading: 13.18-13.20 +/- 0.1, with 1.2 effective runs.
  - Calibration weights change nothing. Survival reading and M92 3-sigma rule: zero weight (0 of 1000 runs within
    3 sigma; best fit 9.6 sigma above the null).
- **Conclusion:**
  - His own pipeline gives 13.2 Gyr with either seed and with 500 or 1000 runs. That is within half his quoted
    uncertainty of 13.39 +/- 0.39 (DM 15.43-15.44 vs 15.49, E 0.080-0.085 vs 0.08).
  - In both seeds the valley bottom is 13.0-13.4 Gyr, and 12.0-12.8 sits 0.2-0.7 noise units higher. So the
    0.4 Gyr between this and our pipeline with his pieces (12.8, `../m15_yingiso_err99`) looks partly systematic, if
    small. It is not only noise in the minimum.
  - His +/- 0.39 Gyr still cannot be reproduced: every fit is about 9.5 sigma outside his null, and his weighting
    concentrates on one isochrone.
