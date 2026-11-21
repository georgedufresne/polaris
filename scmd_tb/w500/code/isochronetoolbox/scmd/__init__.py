"""Synthetic CMDs for NGC 6791 and their comparison with the observed CMD.

``model``     isochrone -> photometry -> one synthetic CMD
``stats``     KDE Jensen-Shannon distance and the 2D KS statistic
``pipeline``  the (age, DM, E(B-V)) grid, units of work for SLURM, collection
``plot``      observed CMD next to an sCMD

Run it with ``python -m isochronetoolbox.scmd`` (see ``__main__``).
"""

from .model import BCTable, Quota, Scatter, apparent, make_scmd, read_isochrones, select_age
from .pipeline import Context, Setup, collect, process_age, regenerate, seed, worker
from .stats import KDE, js_adaptive, js_fixed, ks_2d, rank2d
