#!/usr/bin/env python3
"""Main-sequence core width and KS at the v2 best cell for M. Ying's error model
(MakeTestCMD: per-band binned errors, drawn independently) against ours (paired
errors of the nearest artificial stars, smoothed) with and without the PSF scatter."""
import json, sys
from dataclasses import replace
import numpy as np
B = "/dartfs-hpc/scratch/f007058/scmd_tb"
sys.path.insert(0, f"{B}/tuc_grid_v2/code")
from isochronetoolbox.scmd import Setup, context
from isochronetoolbox.scmd.pipeline import _window
from isochronetoolbox.scmd.stats import ks_2d
S0 = Setup(**{k: v for k, v in json.load(open(f"{B}/tuc_grid_v2/setup.json")).items() if k in Setup.__dataclass_fields__})
CELL = (10005, 11.1, 13.44, 0.030)
BINS = [(15.5, 16.4), (17.6, 18.0), (18.0, 18.4), (18.4, 18.8), (18.8, 19.0)]
mad = lambda x: 1.4826 * np.median(np.abs(x - np.median(x)))


def width(p, lo, hi, step=0.05):
    s = p[(p[:, 1] > lo) & (p[:, 1] < hi)]; r = []
    for a in np.arange(lo, hi, step):
        q = s[(s[:, 1] >= a) & (s[:, 1] < a + step), 0]
        if len(q) > 20: r.append(q - np.median(q))
    return mad(np.concatenate(r))


models = {
    "Martin (independent per-band bins), no PSF": dict(ast_errors="independent", extra_sigma_table=()),
    "ours (paired nearest ASTs), no PSF": dict(extra_sigma_table=()),
    "ours + PSF scatter (v2)": dict(),
}
data = {}
for zp in (True, False):
    c = context(replace(S0, zp_correct=zp, extra_sigma_table=()))
    data[zp] = (c.obs3, c.obs3_cdf)
print("core colour width (1.4826 MAD about the ridge), mmag, per F606W bin")
print(f"{'':46s}" + "".join(f"{lo:.1f}-{hi:.1f}".rjust(11) for lo, hi in BINS) + "   KS vs corr / uncorr data")
for zp in (True, False):
    print(f"{'data, ' + ('zero-point corrected' if zp else 'uncorrected'):46s}"
          + "".join(f"{width(data[zp][0], lo, hi)*1e3:11.1f}" for lo, hi in BINS))
for name, kw in models.items():
    ctx = context(replace(S0, **kw))
    _, make = ctx.model(CELL[0], CELL[1])
    ps = [(lambda p: p[_window(p, S0.ks_window)])(ctx.points(make(np.random.RandomState(700 + k), CELL[2], CELL[3]))) for k in range(4)]
    w = [np.mean([width(p, lo, hi) for p in ps]) for lo, hi in BINS]
    ksc = np.mean([ks_2d(data[True][0], p, data[True][1])[0] for p in ps])
    ksu = np.mean([ks_2d(data[False][0], p, data[False][1])[0] for p in ps])
    print(f"{name:46s}" + "".join(f"{x*1e3:11.1f}" for x in w) + f"   {ksc:.4f} / {ksu:.4f}")
