# polaris

Code written on the Dartmouth Polaris/Discovery scratch space (`/dartfs-hpc/scratch/f007058`) for the sCMD age fits
of NGC 6791, M15 and 47 Tuc. It has the same layout as scratch, so `/dartfs-hpc/scratch/f007058/X` in a script means
`X` here. Only code and notes are tracked (Python, shell/Slurm drivers, READMEs, small configs and namelists). Data,
isochrones, plots and logs are not.

- **Model outputs and isochrones** (physics on, Monte Carlo mode) are archived on the Chaboyer lab share in
  `/dartfs-hpc/rc/lab/C/ChaboyerB/george/`. Its `README.md` says what each archive is and how to restore it.
- **Refresh** this copy with `./sync_from_scratch.sh`, then commit. The script never deletes anything.
- **Archiving scripts** are in `archive/`.

## Where things are

| Topic | Folders |
|---|---|
| sCMD toolbox and all M15 / 47 Tuc fit experiments | `scmd_tb/` (each experiment folder has a README; `scmd_tb/*/code/` is the toolbox snapshot that experiment used) |
| M15 current result (DSEP4 with Monte Carlo fixed) | `scmd_tb/m15_final_mc/` |
| M15 DSEP4 runs (namelists with `LMonteCarlo` fixed, batch scripts) | `M15/dsep4/` |
| M15 vs Ying et al. 2025 replica and tests | `scmd_tb/m15_martin_replica/`, `scmd_tb/m15_*` |
| 47 Tuc DSEP4 rerun with Monte Carlo on | `tuc_dsep4_batch/`, `scmd_tb/tuc_d4v2_build/` |
| 47 Tuc DSEP3 vs DSEP4 figures | `scmd_tb/tuc_iso_compare/` |
| NGC 6791 sCMD legs (OG, DSEP3 physics on, emulator, W v1.4) | `scmd_noco_pipeline/`, `scmd_emul_pipeline/`, `scmd_W14_pipeline/`, `scmd_3way_compare.py` |
| NGC 6791 W emulator isochrones | `emul_iso/`, `emul_W/`, `emul_W14/` |
| NGC 6791 DSEP4 tests and audits | `dsep4_batch/`, `audit_ngc6791/` |

Not included: other people's code (Ying's files, the CONF1DENCE-W snapshot, the DSEP4 distribution) and
`isochrone-toolbox`, which has its own repository.
