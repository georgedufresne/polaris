# pyeep — a portable Python port of MESA `make_eep` and `make_iso`

Self-contained, dependency-light re-implementation of the **`make_eep`** and
**`make_iso`** programs from the Fortran [`iso`](https://github.com/dotbot2000/iso)
project. `make_eep` turns [MESA](http://mesa.sourceforge.net) stellar-evolution
*history* files into **EEP** (Equivalent Evolutionary Point) tracks — the uniform
basis of Dotter (2016) §2 — and `make_iso` interpolates within those tracks to
build **isochrones** (§3.2).

> **Method credit.** The EEP method implemented here is due to Aaron Dotter —
> Dotter, A. 2016, *MIST 0: Methods for the Construction of Stellar Isochrones*,
> **ApJS, 222, 8** ([doi](https://doi.org/10.3847/0067-0049/222/1/8),
> [arXiv:1601.05144](https://arxiv.org/abs/1601.05144)).
> **If you use this module in published work, cite Dotter (2016).**
> Full attribution in [CREDITS.md](CREDITS.md).

This module is designed to be **dropped onto a new cluster and run**: no build
step, no Fortran toolchain, no MESA installation. Clone, `pip install numpy`,
edit one input file, go.

`make_eep` and `make_iso` are ported. `make_track`, `make_cmd`, the
colour/bolometric-correction code and the `*_interp_*` programs are **not**
included.

---

## Contents

```
pyeep/                         (inside the isochrone-toolbox repository)
├── README.md                  this file
├── requirements.txt
├── make_eep.py                main driver (port of src/make_eep.f90)
├── make_iso.py                isochrone driver (port of src/make_iso.f90)
├── run_history.py             convenience runner for raw MESA history.data
├── reformat_mesa.py           MCDSEP/Dartmouth -> MESA layout converter
├── pyeep/                     the library
│   ├── eep.py                 primary + secondary EEP finding
│   ├── iso.py                 isochrone construction (EEP -> isochrone)
│   ├── support.py             track I/O, columns, distance metric
│   ├── phases.py              FSPS phase assignment
│   ├── coord.py               fractional EEP coordinate for every model of a track
│   ├── reformat.py            one-value-per-line -> one-model-per-row
│   └── parallel.py            multi-core execution helpers
├── run/                       ← YOUR WORKING DIRECTORY (see below)
│   ├── input.nml              &eep_controls and &iso_controls namelists
│   ├── input.format           column widths (parity with Fortran)
│   ├── input.eep              EEP interval specification
│   ├── input_dsee.nml         input.nml with the optional EEP rules switched on (DSEE tags v1.2, see below)
│   ├── input_dsee_v1.1.nml    the earlier DSEE tags (v1.1)
│   ├── input_dsee.eep         EEP intervals with the extra RGBBase primary
│   ├── input_dsee_ladder.nml  the Xc/HeL ladder (DSEE tags v1.3, see below)
│   ├── input_dsee_ladder.eep  EEP intervals for the ladder's 19 primaries
│   ├── my_history_columns_basic.list
│   ├── my_history_columns_full.list
│   ├── input.example          ← template, edit this
│   └── input.example.reference  a real known-good config, for comparison
├── tests/
│   ├── run_tests.py           runs the whole self-check suite (no MESA data needed)
│   ├── test_make_eep.py       EEPs from a synthetic history file
│   ├── test_make_iso.py       isochrone kernels + synthetic EEP grid
│   ├── test_pipeline.py       history files -> EEPs -> isochrones, end to end
│   ├── test_golden.py         EEPs and isochrones must match stored reference values
│   ├── data/golden.json       those reference values
│   ├── test_reformat.py       MCDSEP converter
│   ├── test_parallel.py       multi-core helpers
│   ├── test_dsee_options.py   the optional EEP rules of input_dsee.nml
│   ├── test_coord.py          per-model EEP coordinate
│   ├── compare_with_fortran.py      optional: make_eep vs the Fortran binary
│   └── compare_iso_with_fortran.py  optional: make_iso vs the Fortran binary
└── CREDITS.md                 attribution: Dotter (2016), MESA, upstream iso
```

## Requirements

- **Python 3.8+**
- **NumPy**
- `chardet` *(optional)* — better encoding detection in `reformat_mesa.py`;
  without it the reformatter falls back to UTF-8 then latin-1.

```bash
pip install -r requirements.txt
```

No compiler, no MESA, no MESA SDK.

---

## Quick start on a new cluster

```bash
# 1. get the code
git clone https://github.com/georgedufresne/isochrone-toolbox.git
cd isochrone-toolbox/pyeep

# 2. dependencies
pip install --user numpy

# 3. verify the install (uses synthetic data — no MESA files needed)
python3 tests/run_tests.py

# 4. point the input file at your tracks
cd run
cp input.example my_run.example
$EDITOR my_run.example        # set the three directory paths + track list

# 5. run: history files -> EEP tracks -> isochrones
../make_eep.py my_run.example -j 0
../make_iso.py my_run.example
```

Step 3 should end with `28/28 passed`. If it does, the port is working;
anything after that is configuration.

---

## ⚠ The one rule: run from a directory that has the control files

`make_eep.py` reads three auxiliary files **from the current working
directory** — not from the script's location, and there is no command-line flag
to point elsewhere. This is inherited from the Fortran driver and is preserved
deliberately.

| file | how it is found |
|---|---|
| `input.nml` | literally `./input.nml` |
| `input.format` | literally `./input.format` |
| `input.eep` | `./` + whatever `eep_input_file` says in the namelist |
| `my_history_columns_*.list` | `./` + whatever the input file's line 10 says |

That is why the control files live in **`run/`**. The intended invocation is:

```bash
cd run
../make_eep.py my_run.example
```

`make_eep.py` inserts its own directory on `sys.path`, so `pyeep/` is still
imported correctly when you call it by relative path from `run/`.

> **Failure mode to know about.** A missing `input.eep` warns loudly and falls
> back to default EEP intervals. A missing **`input.nml` is silent** — the
> namelist reader returns an empty dict, and `do_phases` then defaults to
> `True` even though the shipped `input.nml` sets `.false.`. You would get
> output, and it would be subtly wrong. If results look unexpected, check you
> are in `run/` first.

---

## The input file

`make_eep.py` takes the same two-part input file as the Fortran program. The
format is **positional** — line order matters, comment lines are counted, so
add and remove nothing:

```
#version string, max 8 characters
pyeep                                     <- max 8 chars
#initial Y, initial Z, [Fe/H], [alpha/Fe], v/vcrit (space separated)
   0.2325 1.44500e-04   -1.88        0.00     0.0
#data directories: 1) history files, 2) eeps, 3) isochrones
/ABSOLUTE/PATH/TO/history_files           <- input tracks live here
/ABSOLUTE/PATH/TO/eeps                    <- .eep output written here
/ABSOLUTE/PATH/TO/isos                    <- isochrones written here by make_iso.py
#read history_columns
my_history_columns_basic.list             <- resolved relative to CWD
#specify tracks
2                                         <- must match the count below
REPLACE_WITH_TRACK_FILENAME_1
REPLACE_WITH_TRACK_FILENAME_2
#specify isochrones
...                                       <- read by make_iso.py (see below)
```

The three directory paths are **absolute** and are the main thing to change
when moving clusters. `input.example.reference` is a real working config from
the Dartmouth cluster (35 tracks, [Fe/H] = −1.88) — useful as a
known-good shape to diff against.

For every listed track the program reads the history file from the *history*
directory, finds primary and secondary EEPs, optionally assigns FSPS phases,
and writes `<name>.eep` into the *eeps* directory.

---

## Running in parallel

Tracks are fully independent — separate input file, separate output file, no
shared mutable state — so `make_eep.py`, `run_history.py` and
`reformat_mesa.py` all accept `-j/--jobs`. (`make_iso.py` takes `-j` too, to
run several grids at once; see below.)

```bash
../make_eep.py my_run.example -j 0     # auto: all cores, capped by free memory
../make_eep.py my_run.example -j 8     # exactly 8 workers
../make_eep.py my_run.example -j -2    # all cores except 2
../make_eep.py my_run.example          # default: serial (unchanged behaviour)
```

The default is **1 (serial)**, so existing invocations behave exactly as
before; with `-j 1` the work runs inline with no pool and no pickling.

Output is **byte-identical regardless of worker count**, and the log stays in
input order: workers buffer their output and the parent replays it in sequence,
so runs remain reproducible and diffable.

**`-j 0` is memory-aware.** Each worker peaks near 190 MiB on ~90 MB history
files, so auto-detect divides available RAM by a 300 MiB allowance and reports
when memory rather than CPU is the limit:

```
note: limiting to 29 workers (9865 MiB available, ~300 MiB per worker); override with -j N
```

On a new cluster this matters: if you are inside a Slurm allocation, `-j 0`
sees the machine's cores, not your cgroup's. Prefer an explicit
`-j $SLURM_CPUS_PER_TASK`.

### Why processes, not threads

The parsing that dominates runtime is pure-Python `str.split`/`float()`, which
holds the GIL, so threads do not overlap — they are *slower than serial* once
contention and overhead are counted:

| backend | wall (35 tracks) |
|---|---|
| serial | 19.8 s |
| `--backend threads -j 16` | 31.1 s (**1.6× slower than serial**) |
| `--backend processes -j 16` | 2.1 s (9.4× faster) |

`--backend threads` exists only so the above can be reproduced. Do not use it
for speed.

---

## Performance (`make_eep`)

Measured on a 35-track grid (80 k–100 k models per track, ~2.9 M model-rows),
16 physical / 32 logical cores:

| | wall clock | speedup |
|---|---|---|
| Fortran `make_eep` (gfortran 12.1, `-O3`) | 46.6 s | 1.0× |
| this port, serial (default) | 18.2 s | 2.6× |
| this port, `-j 0` (auto, 29 procs) | 2.3 s | **20×** |

Scaling by worker count:

| `-j` | 1 | 2 | 4 | 8 | 12 | 16 | 24 | 32 |
|---|---|---|---|---|---|---|---|---|
| wall | 19.8 s | 9.4 s | 5.1 s | 3.0 s | 2.1 s | 2.1 s | 1.9 s | 2.0 s |

Near-linear to ~8 workers, saturating around the 16 physical cores.

Both implementations are dominated by history-file parsing — 93 % of the serial
Python runtime is in `read_history_file`. The Fortran is slower there because it
reads each file twice (once to count rows, once to load) and performs one
internal formatted `read` per column per row (~64 M `read` statements for this
grid). The port tokenises each line once and converts in bulk. The EEP-finding
and interpolation maths — the part that actually differs algorithmically — is a
negligible share of both runtimes.

---

## Making isochrones (`make_iso.py`)

`make_iso.py` reads the `.eep` files written by `make_eep.py` and builds
isochrones by the method of Dotter (2016) §3.2: for every EEP, select the
tracks that define it, interpolate in the age–initial-mass relation to get the
initial mass at the requested age, then interpolate every other column in
initial mass. It takes **the same input file** as `make_eep.py` — the part
after `#specify isochrones` that `make_eep` ignores:

```
#specify isochrones
isochrones_example.txt        <- output name, written into the isochrones directory
min_max                       <- 'min_max' (evenly spaced) or 'list'
linear                        <- 'linear' (years) or 'log10'
81                            <- number of isochrones
4.0e9                         <- min age      (with 'list': one age per line)
12.0e9                        <- max age
double                        <- optional: allow two points per EEP (§3.4)
```

As with `make_eep.py`, run it from `run/`: it reads the **`&iso_controls`**
group of `./input.nml` and the history column list named on line 10 from the
current directory. Unlike `make_eep.py`, a missing `input.nml` is an **error**,
not a silent fallback to defaults. `--nml PATH` points elsewhere.

`&iso_controls` (defaults are the Fortran's; the shipped `run/input.nml`
overrides several):

| parameter | default | meaning |
|---|---|---|
| `do_linear_interpolation` | `.false.` | linear instead of Steffen (1990) monotonic cubic interpolation, in both age→mass and mass→columns |
| `do_smooth` | `.true.` | inverse-distance smoothing of each EEP's age–mass relation before age→mass |
| `do_PAV` | `.true.` | enforce monotonic initial mass along the isochrone (pool adjacent violators, §3.2) |
| `log_age_delta` | `1.0` | only tracks within ±this (dex) of the target age join the interpolation at an EEP |
| `very_low_mass_limit` | `0.5` | tracks at or below this mass are kept even if they lack later primary EEPs, and smoothing is skipped at any EEP that uses one. Set it in `&iso_controls`; the `&eep_controls` value is not read here |
| `set_max_eep_number`, `new_max_eep_number` | `.false.`, `-1` | stop isochrones at this EEP |
| `iso_debug` | `.false.` | print each track's mass and primary EEPs as it is read |

### Many grids at once

Give several input files — one per grid — and `-j` runs them in parallel.
Output is byte-identical to running them one at a time:

```bash
cd run
../make_iso.py /path/to/input_*.example -j 0
```

Measured on 64 of the NGC 6791 grids (39 tracks, 81 ages each), 16 cores:

| | wall clock |
|---|---|
| Fortran `make_iso`, one grid after another | 26.7 s |
| this port, one grid | 0.27 s (Fortran: 0.39 s) |
| this port, `-j 0`, all 64 | **1.5 s (18×)** |

The isochrone arithmetic itself is ~0.6 ms per age: every EEP of an age is
handled in one vectorised pass, and each cubic is evaluated from only the four
tracks around the target rather than rebuilt over the whole grid per column,
as the Fortran does. What remains is mostly writing the text output.

### From Python

```python
# with pyeep/ on sys.path, from a directory holding input.nml (e.g. run/)
import make_iso
isos = make_iso.make_isochrones('my_run.example', write=False)
iso = isos[0]
iso.eep, iso.column('initial_mass'), iso.column('log_L')
```

or, with tracks already in hand, `pyeep.iso.IsochroneGrid(tracks, col_names,
IsoControls(...)).isochrone(age)` — the grid is read once and reused for every
age.

### Where the port deliberately differs from the Fortran

Everything else reproduces the Fortran bit-for-bit, including its quirks, so
that results stay comparable with isochrones already made with it. Two
behaviours are worth knowing about:

- **Repeated interpolation x-values — the port drops the EEP.** When several
  tracks share exactly the same age at an EEP (e.g. very-low-mass tracks that
  all stop at 20 Gyr) MESA's Steffen routine rejects the input and returns
  without computing slopes. The Fortran ignores that error code and evaluates
  coefficients **left in memory by the previous EEP's interpolant**. This was
  confirmed by rebuilding that previous interpolant, which reproduces the
  Fortran's value to all 17 digits. The port drops the point and prints
  `EEP N dropped: repeated x values ...`. It happens only in cubic (non-linear)
  mode; on the NGC 6791 grids it is EEP 353, one point per isochrone. With
  PAV on, the Fortran's stale mass also drags the neighbouring EEPs down when
  pooled (for example, about ten EEPs pulled from 0.74–0.79 to 0.733 M☉),
  so the port's isochrone differs over that short run too. Feeding the Fortran's
  stale value into the port's PAV reproduces the Fortran rows exactly, so
  nothing else is involved.
- **Kept as in the Fortran, so you know:** the "drop non-monotonic points" filter
  (`do l=k-1,1` in `do_isochrone_for_age`) only ever compares the two
  lowest-mass tracks, because that loop runs only when `k = 2`. The mass-order
  check on reading tracks compares each track only with the first one. The port
  reproduces both. The later `masses not monotonic` check still catches an
  out-of-order grid.

---

## Other entry points

### `run_history.py` — skip the input file

For pointing the pipeline at a **real MESA `history.data`** without hand-writing
the two-part input file:

```bash
./run_history.py /path/to/LOGS/history.data
./run_history.py ~/MESA_runs/*/LOGS -o output          # many runs at once
```

It derives the column list from each file's own header row, reads
`&eep_controls`/EEP intervals from `input.nml`/`input.eep` when present, and
writes into `output/` by default. Unlike `make_eep.py` it **does** take
`--nml` and `--eep-input` flags, so it can be run from anywhere:

```bash
./run_history.py /path/to/LOGS/history.data \
    --nml run/input.nml --eep-input run/input.eep \
    --version mist08 --initial-Y 0.2703 --initial-Z 1.42857e-02
```

The derived column list is written per track as `history_columns.<name>.list` —
qualified by track name because columns may differ between runs, and a single
shared file would be a race under `-j`.

### `reformat_mesa.py` — MCDSEP/Dartmouth input

Some generators write history data with **one value per line**: after a
six-line header the file is a flat stream of numbers. `make_eep` expects one
model per row, so transpose first:

```bash
./reformat_mesa.py m0700.feh031.20001              # one file, in place
./reformat_mesa.py /path/to/TEMP20001 -o /path/to/tracks_mesa
```

Inputs may be files or directories (scanned non-recursively). Files already in
MESA layout are detected and skipped, so re-running over a mixed directory is
safe. A model cut off mid-write is dropped rather than emitted as a short row.

Also available programmatically:

```python
from pyeep.reformat import reformat_mesa_file
nmodels = reformat_mesa_file('m0700.feh031.20001', 'track_M0.70_mesa')
```

Full workflow:

```
MCDSEP tracks --reformat_mesa.py--> MESA history --make_eep.py--> .eep --make_iso.py--> isochrones
```

---

## Optional primary-EEP rules (`run/input_dsee.nml`)

`&eep_controls` takes a few switches that change how primary EEPs are placed.
Their defaults reproduce the Fortran exactly; the golden values and the
Fortran comparisons above are unaffected. `run/input_dsee.nml` turns them all
on. These are the rules used to tag the DSEE/CONF1DENCE training grid:

| switch | default (Fortran) | `input_dsee.nml` | effect |
|---|---|---|---|
| `IAMS_Xc` | `0.35` | `0.30` | central X that defines the IAMS (0.30 is MIST's) |
| `RGBTip_require_He_ignition` | `.false.` | `.true.` | no RGB tip unless log(L_He/L_H) > `He_ignition_log_ratio` after TAMS, so a run that stops before the He flash gets no tip instead of one at its last giant-branch model |
| `He_ignition_log_ratio` | `-3.0` | `-3.0` | threshold for the switch above |
| `RGBTip_take_later` | `.false.` | `.true.` | tip = the *later* of max log L and min log Teff (Fortran: the earlier). Metal-poor 2–4 M☉ stars ignite He fainter than they crossed the Hertzsprung gap, so max L alone puts their tip in the gap |
| `add_RGBBase` | `.false.` | `.true.` | extra primary between TAMS and RGBTip at the subgiant-to-giant knee: along TAMS → tip, with x = fraction of the drop to the coolest Teff and y = fraction of the rise to the tip log L, the row that maximises x − y. Needs one more line in the EEP-interval file (`input_dsee.eep`: TAMS at row 454, RGBBase at 491, RGBTip at 605) |
| `VLM_fallback_min_age` | `1.5e10` | `2.09e10` | a star at or below `very_low_mass_limit` whose main sequence outlives the run gets TAMS at its last row only if the run reached this age (the DSEE grid stops at 21 Gyr). The track's `vlm_fallback` flag records it |
| `merger_fallback` | `.true.` | `.false.` | mergers get TAMS at their last row |
| `Xc_HeL_ladder` | `.false.` | `.false.` | a different primary set entirely; see below. `run/input_dsee_ladder.nml` turns it on |

`input_dsee.nml` also sets `logL_scale = 2.0` (the shipped `input.nml` uses
0.125), which weights log L like log Teff in the distance metric that spaces
the secondary EEPs.

To use it, copy `input_dsee.nml` to `input.nml` and keep `input_dsee.eep` next
to it in your working directory, or pass both to `run_history.py`
(`--nml run/input_dsee.nml --eep-input run/input_dsee.eep`).

With `add_RGBBase` every secondary EEP between TAMS and the tip moves, so EEP
numbers from the two configurations are not comparable. Do not mix their
`.eep` files in one isochrone grid.

### The Xc/HeL ladder (`run/input_dsee_ladder.nml`)

`xc_hel_ladder = .true.` replaces everything above the ZAMS with nineteen
primaries that are rungs of the burning rather than features of the HR
diagram:

| # | primary | condition | `.eep` row |
|---|---|---|---|
| 1 | PreMS | first model | 1 |
| 2 | ZAMS | max log g before X_c has dropped 0.001 | 202 |
| 3–13 | X_c ladder | X_c / X_c(ZAMS) ≤ 0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10, 0.01 | 213, 224, 243, 262, 281, 300, 319, 338, 363, 394, 429 |
| 14 | TAMS | X_c ≤ 1e-12 | 454 |
| 15–19 | He ladder | log(L_He / L) ≥ -10, -8, -6, -4, -2 | 479, 492, 517, 554, 605 |

IAMS, RGBBase and RGBTip are not placed, so the four switches that place them
(`IAMS_Xc`, `RGBTip_require_He_ignition`, `RGBTip_take_later`, `add_RGBBase`)
do nothing while the ladder is on. `do_phases` is forced off: the FSPS phases
are defined by the standard primaries.

Why: the two anchors that carried the post-main-sequence in 1.0–1.2 are
morphological, and both fail on part of the grid. The RGB tip is the later of
max log L and min log Teff, so a track that stops before the He flash has no
tip at all, and a track that stops just past a local maximum gets one in the
wrong place -- on the DSEE grid roughly a quarter of the tracks. RGBBase is a
knee in the HR diagram whose position depends on the whole TAMS → tip segment,
so it moves when the tip does. The ladder's rungs are thresholds on X_c and on
L_He/L, both monotone, both local: a rung is either reached or it is not, and
a truncated track simply has no EEP above the rung it reached, instead of a
spurious one at its last row. The same rung is the same stage of burning at
every mass, which is what an emulator conditioned on EEP number needs.

The intervals are not uniform. The main sequence is sampled lightly where the
track barely moves (0.95 → 0.80, ten to eighteen points per rung) and heavily
as the core runs out (0.20 → 0.01, twenty-four to thirty-four), and the giant
branch is weighted towards the top (fifty points in the last rung), where the
spacing in Teff and log L is largest. TAMS still lands on row 454 and the last
primary on row 605, as in `input_dsee.eep`.

The He rungs read `log_LHe - log_L`, i.e. genuine log10 columns in the MESA
sense. DSEP writes those columns linearly in some output formats; with linear
values every rung fires at the TAMS.

To use it: `--nml run/input_dsee_ladder.nml --eep-input run/input_dsee_ladder.eep`.

### Versions of the DSEE tags

Each set of DSEE/CONF1DENCE training data was tagged with one of these
configurations. Each reproduces the tagger that made the data exactly: same
primary EEPs and same per-model coordinates on 2,922 grid tracks (10.8 M models).

| version | configuration | used for | what changed |
|---|---|---|---|
| 1.0 | `input.nml` + `input.eep` (the Fortran rules) | first EEP grid (emulators A–J, C, Cx) | -- |
| 1.1 | `input_dsee_v1.1.nml` + `input.eep` | retagged grid (emulator R) | RGB tip only after He ignition and at the later of max L / min Teff; IAMS at X_c = 0.30; flagged very-low-mass fallback; no merger fallback |
| 1.2 | `input_dsee.nml` + `input_dsee.eep` | current grid (emulators S, V, W, T) | 1.1 + RGBBase primary (row 491); log L weighted like Teff in the distance metric (`logL_scale` 2.0); very-low-mass fallback only for runs that reached 21 Gyr |
| 1.3 | `input_dsee_ladder.nml` + `input_dsee_ladder.eep` | -- | the Xc/HeL ladder: 19 primaries, no IAMS, RGBBase or RGBTip |

EEP numbers are not comparable between versions: changing a primary or the
distance metric moves every secondary EEP in that interval (1.1 → 1.2 moves all
of rows 454–605; 1.3 changes every row above 202). Compare tracks or isochrones from different versions in
physical quantities (Teff, L, age), never by EEP number.

### Fractional EEP coordinate for every model

`make_eep` writes resampled tracks. For data that keeps every MESA model (e.g.
emulator training sets), `pyeep.coord.eep_coordinate(track, state)` gives each
model its fractional `.eep` row: exactly the primary's row at each primary,
non-decreasing in between, NaN after the last primary. Reading a history column
at integer coordinates reproduces the `.eep` track (tested in
`tests/test_coord.py`).

```python
from pyeep.coord import eep_coordinate
coord, n_primary, primary_rows = eep_coordinate(track, state)   # same controls as make_eep
```

## Differences from the Fortran `make_eep`

(For `make_iso`, see [Where the port deliberately differs](#where-the-port-deliberately-differs-from-the-fortran)
above.) Deliberate, behaviour-preserving simplifications:

- **NumPy instead of MESA libraries** for array math and interpolation.
- **History parsing by whitespace tokenising** rather than fixed-width column
  slicing — a robust equivalent for standard MESA output, yielding identical
  values; requested columns are matched against the column-name header row.
- **No binary (`.bin`) caching.** The Fortran can cache tracks/EEPs in
  unformatted binary; the port always reads/writes text.
- The **cubic-interpolation branch** in `eep_interpolate` is dead code in the
  Fortran (guarded by `linear = .true.`) and is not ported.
- **`libm` `log10` rather than NumPy's `log10`** in `distance_along_track`. The
  two disagree by 1 ULP on ~5 % of inputs; because the distance metric is a
  cumulative sum over the whole track (often >10⁵ steps), those ULPs accumulate
  and shift every interpolated secondary EEP by ~10⁻¹⁰ relative. Using
  `math.log10` (the same libm entry point gfortran calls) reproduces the
  Fortran distances bit-for-bit.

### One intentional behavioural divergence: `initial_mass`

In the Fortran, `track% initial_mass` has **no default initialiser**, and
`read_history_file` calls `set_star_type_from_history` *before* recovering the
mass from the track data. When the MESA header carries no `initial_mass` field
— as with standard `history.data` — the sub-stellar test

```fortran
if( t% initial_mass < 0.1d0 .and. maxval(abs(t% tr(i_Xc,:) - t% tr(i_Xc,1))) < 0.1d0 )
```

reads **uninitialised heap memory**. Since `make_eep.f90` allocates and
deallocates `t` once per track, that memory is normally the *previous* track's
`initial_mass` (garbage, usually 0, for the first track). A track is then
mislabelled `substellar`, which makes `ZAMS` return `argmax(log_center_T)`
instead of `argmax(log_g)` — a completely different EEP.

This only bites tracks whose central hydrogen changes by less than 0.1 over the
whole track. It is demonstrable: a 0.24 M☉ track placed **first** in the input
list gets `ZAMS = 5115` from the Fortran, but `ZAMS = 79` when placed
**second**, behind a 0.28 M☉ track — same file, same controls.

The port resolves `initial_mass` from the track data *before* classifying, so
it is deterministic and independent of input ordering. It always produces
`ZAMS = 79`, which is also what the Fortran produces whenever the preceding
track has mass ≥ 0.1. A 0.24 M☉ star is a hydrogen-burning M dwarf, not a brown
dwarf, so `low-mass` is the intended classification.

**Practical consequence:** `.eep` files for very low-mass tracks generated by
the Fortran may differ from this port's — header `substellar` vs `low-mass`,
and every data row shifted. The port's answer is the correct one.

---

## Validation

### Self-check (ships with the package)

One command runs every test and prints a pass/fail summary (about a second;
covers the EEP side, the isochrone side, and `tests/test_pipeline.py`, which
makes both in sequence — synthetic history files → `make_eep` → `make_iso` —
and checks the recovered masses against a known age–mass relation;
exit status non-zero on failure).

`tests/test_golden.py` pins the actual numbers: EEP tracks and isochrones made
from synthetic inputs are compared with `tests/data/golden.json`, produced by
the version of the code validated against the Fortran binaries. A deliberate
change to the numerics needs `PYEEP_REGEN_GOLDEN=1 python3 tests/test_golden.py`
and a reason in the commit message.

The suite was checked by planting bugs on purpose; each one fails it:

| planted bug | caught by |
|---|---|
| Steffen coefficient off by 0.05 % | `test_make_iso`, `test_golden` |
| PAV switched off in isochrones | `test_make_iso::test_pav_enforces_monotonic_mass` |
| EEP distance metric missing its square root | `test_golden` (EEP and isochrone cases) |
| secondary EEPs spaced off by one | `test_make_eep`, `test_pipeline`, `test_golden` |

Passes on Python 3.8, 3.9 and 3.12 (NumPy 1.24 – 2.5), including on 2 CPUs. GitHub runs the same command on every push
that touches `pyeep/` — see `.github/workflows/pyeep-tests.yml`.

```bash
python3 tests/run_tests.py          # everything
python3 tests/run_tests.py iso      # only tests with 'iso' in the name
```

Or file by file:

```bash
python3 tests/test_make_eep.py     # end-to-end on a synthetic MESA track
python3 tests/test_make_iso.py     # isochrone kernels and a synthetic EEP grid
python3 tests/test_pipeline.py     # history files -> EEPs -> isochrones
python3 tests/test_golden.py       # reference values for EEPs and isochrones
python3 tests/test_reformat.py     # MCDSEP converter round-trip
python3 tests/test_parallel.py     # worker-count, ordering, error propagation
```

No MESA data required — synthetic history files and EEP grids are generated on
the fly.

Optionally point at a real file for an extra regression:

```bash
PYEEP_REAL_HISTORY=/path/to/LOGS/history.data python3 tests/test_make_eep.py
```

### Numeric regression against the Fortran binary

`tests/compare_with_fortran.py` runs both implementations over the same grid
and diffs the results. It needs a compiled Fortran `make_eep`, which you must
build yourself from [upstream `iso`](https://github.com/dotbot2000/iso) and
which in turn needs a full MESA install — so it is optional, and not part of
the self-check. Upstream sources are deliberately not vendored here; see
[CREDITS.md](CREDITS.md).

Results on a 35-track, [Fe/H] = −1.88 grid (0.24–2.40 M☉, MESA r15140):

| | |
|---|---|
| Tracks byte-identical to the Fortran | **34 / 35** |
| Worst relative difference | **2.2 × 10⁻¹⁶** (1 ULP, 3 values in one track) |
| Primary-EEP indices identical | **35 / 35** |
| Rows copied verbatim at primary EEPs | bit-exact |

The one non-matching track is the 0.24 M☉ case above — the expected and
intended difference.

On a second, independent grid — 39 MCDSEP tracks (`feh031`, id 20001) passed
through `reformat_mesa.py` first, exercising converter and port together:

| | |
|---|---|
| Tracks matching the Fortran | **38 / 39** |
| Worst relative difference | **7.2 × 10⁻¹⁶** |

The one exception is the 0.20 M☉ track, the first in the list: the
`initial_mass` case described above (Fortran `substellar`, port `low-mass`).

The converter alone was byte-identical to the original standalone reformatting
script in all 39 cases.

### `make_iso` against the Fortran binary

`tests/compare_iso_with_fortran.py` runs the Fortran `make_iso` and this port
on the same `.eep` grids under six `&iso_controls` configurations, then
compares the output files:

```bash
python3 tests/compare_iso_with_fortran.py --fortran /path/to/make_iso \
    --columns run/my_history_columns_basic.list input_A.example input_B.example
```

Results on 13 randomly chosen NGC 6791 grids (81 ages each):

| configuration | result |
|---|---|
| linear interpolation, no smoothing/PAV (the shipped `input.nml`) | **13/13 byte-identical** |
| `set_max_eep_number` | **13/13 byte-identical** |
| `double` EEP mode | **13/13 byte-identical** |
| Steffen cubic | identical except the dropped EEP 353 (see above) |
| Steffen + smoothing + PAV (Fortran defaults), linear and log10 ages | identical except EEP 353 and, with PAV, its pool |

A further 64 grids at the shipped settings were **64/64 byte-identical** via
`-j 0`. End to end on grid 20001 (history files → `make_eep` → `make_iso`),
the Python chain is byte-identical to the Fortran chain once both use the same
0.20 M☉ EEP track.

The self-check `tests/test_make_iso.py` needs no data: it compares the
vectorised Steffen interpolation, binary search and output formatting with
direct scalar transcriptions of the Fortran/MESA routines on thousands of
random inputs (bit-identical), and recovers the initial masses of a synthetic
grid built from a known age–mass relation.

---

## Notes on real MESA data

Two compatibility points, both handled:

- **String version numbers.** Older MESA wrote an integer `version_number`;
  newer releases write a quoted string such as `"r22.05.1"`. The reader accepts
  both, extracting the leading integer for the `.eep` header.
- **Missing optional columns.** A run truncated early may omit `center_gamma`
  or `surface_h1`. These are handled gracefully — the affected EEP tests and
  phase promotions are skipped — matching the fact that such a run cannot reach
  the evolutionary points that need them.

---

## Relationship to the Fortran original

This module re-implements five files from Dotter's `iso`, plus the
piecewise-monotonic interpolator of MESA's `interp_1d` used by `make_iso`:

| Fortran source (upstream `iso`) | Python module here |
|---|---|
| `src/make_eep.f90` | `make_eep.py` |
| `src/make_iso.f90` | `make_iso.py`, `pyeep/iso.py` |
| `src/eep.f90` | `pyeep/eep.py` |
| `src/phases.f90` | `pyeep/phases.py` |
| `src/iso_eep_support.f90` | `pyeep/support.py`, `pyeep/iso.py` (relevant subset) |
| MESA `interp_1d` (`interp_pm`, `interp_value`) | `pyeep/iso.py` (`steffen_eval`) |

`src/eep_approx.f90` and `src/eep_exact.f90` are not part of `make_eep` (per the
upstream makefile) and were not ported.

`pyeep/reformat.py` + `reformat_mesa.py` and `pyeep/parallel.py` have no Fortran
counterpart — they are additions.

**No upstream source is redistributed in this repository.** The `iso` project
carries no licence file, so its sources remain all-rights-reserved. To build the
Fortran for comparison, clone it yourself:

```bash
git clone https://github.com/dotbot2000/iso    # needs MESA, $MESA_DIR set
```

See [CREDITS.md](CREDITS.md) for the full citation and attribution.
