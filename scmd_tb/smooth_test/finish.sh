#!/bin/bash
B=/dartfs-hpc/scratch/f007058/scmd_tb/smooth_test
PY=/dartfs-hpc/scratch/f007058/envs/fidanka-env/bin/python
export HOME=/dartfs-hpc/scratch/f007058
until [ "$(ls $B/smooth_resample/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 200 ] && [ "$(ls $B/smooth_noresample/bootstrap/part_*.csv 2>/dev/null | wc -l)" -ge 40 ] && [ "$(ls $B/ying_null/parts/part_*.csv 2>/dev/null | wc -l)" -ge 31 ]; do sleep 30; done
sleep 20; scancel -n smoothA -u f007058; scancel -n smoothB -u f007058; scancel -n yingnull -u f007058
for d in smooth_resample smooth_noresample; do cd $B/$d/code && echo "== $d" && $PY -m isochronetoolbox.scmd bootstrap-collect --out $B/$d 2>&1 | grep -v Warn; done
$PY - <<'PY'
import glob, pandas as pd
d = pd.concat([pd.read_csv(f) for f in glob.glob('/dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/ying_null/parts/part_*.csv')])
d.to_csv('/dartfs-hpc/scratch/f007058/scmd_tb/smooth_test/ying_null/ying_null.csv', index=False)
for (c, s), g in d.groupby(['config', 'size']):
    print(f'Ying-style null [{c}, resample {s}]: n {len(g)}, KS {g.ks.mean():.4f} +- {g.ks.std(ddof=0):.4f}, 3 sigma {g.ks.mean() + 3 * g.ks.std(ddof=0):.4f}')
PY
