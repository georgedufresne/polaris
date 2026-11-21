#!/usr/bin/env python3
"""Write Ying et al. (2025)'s M15 isochrones for our 500 runs in the exact block format that M. Ying's MakeTestCMD.f90
reads, as feh227cmd.<run> (the name his TestCMDPAR.sh passes; MakeTestCMD takes the run number from the 5 characters
after "cmd.").

Source: /dartfs-hpc/scratch/f007058/M15/iso_ying2025/<run>/ying_cmd.<run> (his Zenodo HDF5 values, written and
validated exactly by ../m15_yingiso_nozp/convert_ying_isochrones.py).

MakeTestCMD's reader, per age block: one skipped line; a header read with format (1x,I3,F9.6,F8.4,F7.0) giving npts,
mixing length, ALPHAC, age (Myr); one skipped line; npts list-directed rows "EEP mass F606W F606W-F814W"; two skipped
lines. The ying_cmd files have no trailing lines and a header whose fields do not sit in those columns (the age would be
read as 18000), so they cannot be passed to it directly.

Validation: every header is parsed back with the same fixed columns (npts and age must match), and every data row is
identical to the source. Output: /dartfs-hpc/scratch/f007058/M15/iso_ying2025_martinfmt/, logs/write_martin_isochrones.out."""
import os

SRC = "/dartfs-hpc/scratch/f007058/M15/iso_ying2025"
OUT = "/dartfs-hpc/scratch/f007058/M15/iso_ying2025_martinfmt"
E = "/dartfs-hpc/scratch/f007058/scmd_tb/m15_martin_replica"
os.makedirs(OUT, exist_ok=True)


def blocks(path):
    """(npts, mixlen, alphac, age_myr, rows) per age block of a ying_cmd file."""
    lines = open(path).read().splitlines()
    i, out = 0, []
    while i < len(lines):
        assert lines[i].startswith("#NPTS"), (path, i, lines[i])
        f = lines[i + 1][1:].split()
        n, ml, ac, age = int(f[0]), float(f[1]), float(f[2]), float(f[3])
        rows = lines[i + 3:i + 3 + n]
        assert len(rows) == n and not any(r.startswith("#") for r in rows)
        out.append((n, ml, ac, age, rows))
        i += 3 + n
    return out


def fortran_header(s):
    """(1x,I3,F9.6,F8.4,F7.0) with blanks ignored (Fortran's default BLANK='NULL')."""
    def num(t, conv):
        return conv(t.replace(" ", ""))
    return (num(s[1:4], int), num(s[4:13], float), num(s[13:21], float), num(s[21:28], float))


runs = [int(x) for x in open(f"{E}/runs.txt").read().split()]
nblk = 0
for r in runs:
    bl = blocks(f"{SRC}/{r}/ying_cmd.{r}")
    txt = []
    for n, ml, ac, age, rows in bl:
        assert n <= 999
        hdr = f"#{n:3d}{ml:9.6f}{ac:8.4f}{age:6.0f}."
        assert len(hdr) == 28
        txt += [f"#NPTS {n}", hdr, "# EEP Mass F606W F606W-F814W", *rows, "", ""]
    p = f"{OUT}/feh227cmd.{r}"
    open(p, "w").write("\n".join(txt) + "\n")
    # read back as MakeTestCMD does
    L = open(p).read().splitlines()
    i = 0
    for n, ml, ac, age, rows in bl:
        hn, hml, hac, hage = fortran_header(L[i + 1])
        assert hn == n and abs(hage - age) < 0.01 and abs(hml - ml) < 1e-6, (r, L[i + 1])
        assert L[i + 3:i + 3 + n] == rows
        assert L[i + 3 + n] == "" and L[i + 4 + n] == ""
        i += 5 + n
        nblk += 1
    assert i == len(L), (r, i, len(L))
print(f"{len(runs)} files, {nblk} age blocks written to {OUT}; every header parses back (npts, age, mixing length) with "
      f"MakeTestCMD's fixed columns and every row matches the source")
