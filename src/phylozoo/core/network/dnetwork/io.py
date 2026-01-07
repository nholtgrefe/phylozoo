"""
Network I/O module.

This module registers format handlers for DirectedPhyNetwork with FormatRegistry
for use with the IOMixin system.

All eNewick functionality (parser, writer, reader) is in the _enewick module.
"""

from __future__ import annotations

from ....utils.io import FormatRegistry
from .base import DirectedPhyNetwork
from ._enewick import to_enewick, from_enewick

# Register eNewick format handlers with FormatRegistry
FormatRegistry.register(
    DirectedPhyNetwork, 'enewick',
    reader=from_enewick,
    writer=to_enewick,
    extensions=['.enewick', '.eNewick', 'enwk', '.nwk', '.newick'],
    default=True
)
