"""pyeep: a Python port of the MESA ``make_eep`` and ``make_iso`` programs.

Transforms MESA stellar-evolution history files into EEP-based tracks, and
EEP tracks into isochrones.  See :mod:`pyeep.support`, :mod:`pyeep.eep`,
:mod:`pyeep.phases` and :mod:`pyeep.iso`.
"""

from .support import (  # noqa: F401
    Column,
    Track,
    EEPState,
    setup_columns,
    set_eep_interval,
    set_default_eep_interval,
    read_history_file,
    write_track,
)
from .eep import primary_eep, secondary_eep  # noqa: F401
from .phases import set_track_phase  # noqa: F401
from .iso import (  # noqa: F401
    IsoControls,
    Isochrone,
    IsochroneGrid,
    read_eep,
    write_isochrones,
)

__all__ = [
    'Column', 'Track', 'EEPState',
    'setup_columns', 'set_eep_interval', 'set_default_eep_interval',
    'read_history_file', 'write_track',
    'primary_eep', 'secondary_eep', 'set_track_phase',
    'IsoControls', 'Isochrone', 'IsochroneGrid', 'read_eep', 'write_isochrones',
]
