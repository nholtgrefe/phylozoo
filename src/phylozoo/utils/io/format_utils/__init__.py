"""Shared format structure (NEXUS, PHYLIP). Class-specific (de)serialization lives in ``core/*/io``."""

from . import nexus
from . import phylip

# Re-export all public functions from nexus and phylip for convenient access
from .nexus import (
    nexus_header,
    parse_nexus,
    write_block,
    write_taxa_block,
)
from .phylip import (
    parse_phylip_matrix,
    write_phylip_matrix,
)

__all__ = [
    'nexus',
    'phylip',
    # Nexus
    'nexus_header',
    'parse_nexus',
    'write_block',
    'write_taxa_block',
    # PHYLIP
    'parse_phylip_matrix',
    'write_phylip_matrix',
]
