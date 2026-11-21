"""2D fits at 11.4 Gyr (all 341 cells of run 10099): KS^2 = K0^2 + quadratic form in (mu, E)."""
import sys
import numpy as np
sys.path.insert(0, '/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/discretisation')
from common import load, Model, STEP, MU0, E0, AGE0

d = load()
a = d[d.src == 'allcells']
X = a[['v', 'w']].to_numpy(float)
y = a.ks.to_numpy()
print('cells at 11.4 Gyr:', len(a))
print('grid min', y.min(), 'at', a.loc[a.ks.idxmin(), ['mu', 'E']].tolist())
for T in [0.03, 0.04, 0.05, 0.06, 0.07, 0.09]:
    sel = y < T
    print(f'\n== KS < {T}: {sel.sum()} cells, E values {sorted(set(a.E[sel].round(3)))}')
    for kind in ['hyp', 'par', 'cone', 'pfree']:
        m = Model(kind, k=2)
        r = m.fit(X[sel], y[sel], n_starts=16)
        K0, x0, A, p = m.unpack(r.x)
        res = r.fun
        rms = np.sqrt(np.mean(res ** 2))
        dof = sel.sum() - m.npar()
        pred_node = m.predict(r.x, np.array([[0., 0.]]))[0]
        print(f'  {kind:5s} K0={abs(K0):.5f} mu0={MU0 + x0[0]*STEP[1]:.4f} E0={E0 + x0[1]*STEP[2]:.5f} '
              f'p={p if p is None else round(p, 2)} rms={rms:.5f} (dof {dof}) pred@node={pred_node:.5f} '
              f'A(steps)=[[{A[0,0]:.3g},{A[0,1]:.3g}],[{A[1,1]:.3g}]]')
