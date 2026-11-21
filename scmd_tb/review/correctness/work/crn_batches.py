import sys; sys.path.insert(0, ".")
from common import *
S = load("dev")
from isochronetoolbox.scmd import ast as AM
s = setup_fresh(S); ctx = S.context(s)
calls = []
orig = AM.ArtificialStars.radii
def spy(self, rs, n, mag=None):
    calls.append(n); return orig(self, rs, n, mag)
AM.ArtificialStars.radii = spy
_, make = ctx.model(10099, 11.4)
for dm in s.dm[[20, 21, 22, 23, 24]]:
    for e in s.ebv[[2, 3, 4]]:
        calls.clear()
        cat = make(np.random.RandomState(s.crn_seed), dm, e)
        print(f"DM {dm:.2f} E {e:.3f}: batch sizes {calls}")
