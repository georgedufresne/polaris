# iso_truncation: truncated ages in the W isochrone sets (2026-11-22)

**Why.** After the rerun, run 20374 of the W v1.4 NGC 6791 leg (`../../2d_ks_W14`) still had 80 of 81 ages. The cause
is its 7.2 Gyr isochrone in `../../Isochrones_W14/20374`, which has 926 points against about 2,320 at 7.1 and 7.3 Gyr.
Every one of the 308 sCMDs for that age failed the KDE step with "Too few stars". The original run had the same
gap, so this is a defect in the isochrone, not a pipeline failure.

The Z-matched isochrone (`../../Isochrones_W14feh/20374`) is complete at 7.2 Gyr (2,318 points).

A partly truncated curve could pass the KDE step and still bias a fit. This scan checks every run of both sets.

**What.** `scan.sh <iso_dir> <label>` (Slurm array, `--account=free`, 16 tasks per set) counts the isochrone points per
age in each run's `isos/isochrones_example.txt`. It flags any age with fewer than 90 % of the mean count of its two
neighbouring ages.

**Output.** `output/<label>_all.tsv` has one line per run: run, number of ages, and either `ok` or the flagged ages as
`age:count/neighbour_mean`. The per-task files are `output/<label>_<task>.tsv`.

**Jobs.** See `jobids.txt`.
