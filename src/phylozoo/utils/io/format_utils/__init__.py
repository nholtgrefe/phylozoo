"""Shared format structure (NEXUS, PHYLIP, DOT). Class-specific (de)serialization lives in ``core/*/io``."""

from . import nexus
from . import phylip
from . import dot

# Re-export all public functions from nexus, phylip and dot for convenient access
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
from .dot import (
    convert_node_id,
    dot_edge_line,
    dot_node_line,
    escape_dot_string,
    format_dot_attributes,
    parse_dot_attributes,
    parse_dot_document,
)

__all__ = [
    "nexus",
    "phylip",
    "dot",
    # Nexus
    "nexus_header",
    "parse_nexus",
    "write_block",
    "write_taxa_block",
    # PHYLIP
    "parse_phylip_matrix",
    "write_phylip_matrix",
    # DOT
    "convert_node_id",
    "dot_edge_line",
    "dot_node_line",
    "escape_dot_string",
    "format_dot_attributes",
    "parse_dot_attributes",
    "parse_dot_document",
]
