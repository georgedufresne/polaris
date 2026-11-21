import sys; sys.path.insert(0, ".")
from common import *
S = load("dev")
which = sys.argv[1]
s = setup_fresh(S, ast_errors=which); ctx = S.context(s)
_, make = ctx.model(10099, 11.4)
ks = [ctx.stage3(make(np.random.RandomState(500000 + j), 13.38, 0.015))["ks_statistic"] for j in range(12)]
print(which, "data-vs-model KS: mean %.5f  sd %.5f  (n=%d)" % (np.mean(ks), np.std(ks, ddof=1), len(ks)), np.round(ks, 4))
