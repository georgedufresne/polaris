import sys; sys.path.insert(0, ".")
from common import *
S = load("dev")
from isochronetoolbox.scmd import ast as AM, ridge as RG
from isochronetoolbox.scmd.ast import centre, F32
radial = sys.argv[1]
s = setup_fresh(S, radial=radial); ctx = S.context(s)
cap = {}
o_radii, o_draw, o_trace = AM.ArtificialStars.radii, AM.NearestErrors.draw, RG.RidgeCut.trace
def sp_radii(self, rs, n, mag=None):
    r = o_radii(self, rs, n, mag); cap.setdefault("rad", []).append(r); cap.setdefault("v", []).append(np.array(mag)); return r
def sp_draw(self, rs, radius, v, i):
    out = o_draw(self, rs, radius, v, i); cap.setdefault("err", []).append(out); return out
def sp_trace(self, c, m):
    rd = o_trace(self, c, m); cap["ridge"] = rd; return rd
AM.ArtificialStars.radii, AM.NearestErrors.draw, RG.RidgeCut.trace = sp_radii, sp_draw, sp_trace
_, make = ctx.model(10099, 11.4)
res = []
for seed in range(3):
    cap.clear()
    cat = make(np.random.RandomState(500000 + seed), 13.38, 0.015)
    rad = np.concatenate(cap["rad"]); v = np.concatenate(cap["v"])
    ok = np.concatenate([e[0] for e in cap["err"]]); ev = np.concatenate([e[1] for e in cap["err"]])
    # observed magnitudes as make_scmd_ast forms them (sign -1), I from the stored draws
    ei = np.concatenate([e[2] for e in cap["err"]])
    # true I is needed: recover from captured inputs of draw via v_inst/i_inst -- re-derive with the spy below
    res.append((rad, v, ok, ev, ei))
AM.NearestErrors.draw = o_draw
