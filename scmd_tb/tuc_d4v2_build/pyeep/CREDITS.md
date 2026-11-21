# Credits and provenance

## Method and original implementation

The EEP (Equivalent Evolutionary Point) method implemented here — the primary
and secondary EEP definitions, the distance metric along a track, and the
interpolation scheme used to place secondary EEPs — is due to **Aaron Dotter**.

> Dotter, A. 2016, *MESA Isochrones and Stellar Tracks (MIST) 0: Methods for
> the Construction of Stellar Isochrones*, **ApJS, 222, 8**
> doi:[10.3847/0067-0049/222/1/8](https://doi.org/10.3847/0067-0049/222/1/8) ·
> arXiv:[1601.05144](https://arxiv.org/abs/1601.05144)

```bibtex
@ARTICLE{Dotter2016,
   author = {{Dotter}, Aaron},
    title = "{MESA Isochrones and Stellar Tracks (MIST) 0: Methods for the
              Construction of Stellar Isochrones}",
  journal = {\apjs},
     year = 2016,
    month = jan,
   volume = {222},
   number = {1},
      eid = {8},
    pages = {8},
      doi = {10.3847/0067-0049/222/1/8},
   eprint = {1601.05144},
   adsurl = {https://ui.adsabs.harvard.edu/abs/2016ApJS..222....8D},
}
```

**If you use this module in published work, cite Dotter (2016).** The method is
his; this module is only a re-implementation of it.

The original Fortran implementation is the `iso` package:
<https://github.com/dotbot2000/iso>

Isochrone construction (`make_iso`) interpolates with the piecewise-monotonic
cubic of **Steffen (1990)**, as implemented in MESA's `interp_1d` module; the
port re-implements that algorithm (no MESA code is copied):

> Steffen, M. 1990, *A simple method for monotonic interpolation in one
> dimension*, **A&A, 239, 443**

## What this module is

`pyeep` is an independent Python re-implementation of the `make_eep` and
`make_iso` programs from `iso`, written to reproduce its numerical output while running faster and
without a MESA/Fortran toolchain. It was validated against the original binary:
34 of 35 tracks byte-identical, worst relative difference 2.2 × 10⁻¹⁶; the
isochrones are byte-identical to the Fortran `make_iso` at the shipped settings
(see the module README).

**No source code from `iso` is redistributed in this repository.** The upstream
project carries no licence file, so its sources are all-rights-reserved and are
deliberately not vendored here. To run the Fortran regression comparison
(`tests/compare_with_fortran.py`) you must obtain and build `iso` yourself:

```bash
git clone https://github.com/dotbot2000/iso
# requires a full MESA installation with $MESA_DIR set
```

## MESA

Tracks consumed by this module are produced by **MESA** (Modules for
Experiments in Stellar Astrophysics), Paxton et al. 2011, 2013, 2015, 2018,
2019. See <https://docs.mesastar.org> for the current citation set.

## This re-implementation

Python ports of `make_eep` and `make_iso`, parallel execution layer, and the MCDSEP/Dartmouth input converter
by **George Dufresne** (Dartmouth College), MIT licensed — see `LICENSE` in this directory.
