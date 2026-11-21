"""47 Tuc mock sCMDs: full DM x E(B-V) grid at one age, best fits plotted against the observed CMD."""
import sys, time
import pandas as pd
sys.path.insert(0, '/dartfs-hpc/rc/home/8/f007058/isochrone-toolbox')
from isochronetoolbox.scmd import preset, context, process_age, regenerate
from isochronetoolbox.scmd.plot import cmd_pair

OUT = '/dartfs-hpc/scratch/f007058/scmd_tb/tuc_mock'
jobs = [(int(a), float(b), {}, 'pyeep + MIST ACS BCs') for a, b in (x.split(':') for x in sys.argv[1:])]
rows = []
for run, age, ov, what in jobs:
    s = preset('47tuc', **ov); ctx = context(s); t = time.time()
    js, ks = process_age(ctx, run, age)
    ks = sorted(ks, key=lambda r: r['ks_statistic'])
    print(f'run {run} age {age}: {len(js)} cells scored in stage 2, {len(ks)} through stage 3, {time.time() - t:.0f} s', flush=True)
    for k, r in enumerate(ks[:3], 1):
        cat = regenerate(ctx, run, age, r['distance_modulus'], r['extinction_ebv'])
        title = (f'47 Tuc mock  run {run}  {age:.1f} Gyr ({what})\nmu {r["distance_modulus"]:.2f}  E(V-I) {r["extinction_ebv"]:.3f}'
                 f'   KS {r["ks_statistic"]:.4f}   JS {r["js_divergence"]:.4f}   stage-2 rank {r["kde_rank"]}')
        path = f'{OUT}/mock_run{run}_age{age:.1f}_best{k}.png'
        (c0, c1), (m0, m1) = s.ks_window
        cmd_pair(ctx.obs, ctx.points(cat), title, path, window=s.ks_window, labels=('F606W - F814W', 'F606W'),
                 xlim=(0.38, 0.86), ylim=(m1 + 0.15, m0 - 0.15))
        rows.append(dict(run=run, age=age, isochrones=what, rank=k, **{c: r[c] for c in
                    ('distance_modulus', 'extinction_ebv', 'ks_statistic', 'js_divergence', 'kde_js', 'kde_rank', 'n_observed', 'n_simulated')}))
        print('   ', path, flush=True)
    pd.DataFrame(js).to_csv(f'{OUT}/kde_js_run{run}_age{age:.1f}.csv', index=False)
pd.DataFrame(rows).to_csv(f'{OUT}/mock_summary.csv', index=False)
print(pd.DataFrame(rows).to_string(index=False))
