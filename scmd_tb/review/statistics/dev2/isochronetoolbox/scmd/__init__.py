"""Synthetic CMDs for NGC 6791 and their comparison with the observed CMD.

``model``     isochrone -> photometry -> one synthetic CMD (Gaia, regional scatter)
``ast``       the same with artificial-star errors and completeness (HST, 47 Tuc)
``ridge``     ridge line through a CMD (G. Dufresne's NGC 6791 method) and the cut around it
``stats``     KDE Jensen-Shannon distance and the 2D KS statistic
``pipeline``  the (age, DM, E(B-V)) grid, units of work for SLURM, collection
``plot``      observed CMD next to an sCMD

Run it with ``python -m isochronetoolbox.scmd`` (see ``__main__``).
"""

from .model import BCTable, Quota, Scatter, apparent, make_scmd, read_dsep, read_isochrones, select_age
from .pipeline import (PRESETS, ASTContext, Context, Setup, collect, context, preset, process_age,
                       regenerate, seed, worker)
from .ridge import Ridge, RidgeCut, Tracer
from .stats import KDE, js_adaptive, js_fixed, ks_2d, rank2d
