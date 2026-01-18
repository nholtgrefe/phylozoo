"""
Visualization module for PhyloZoo (viz).

This module provides a modular visualization system for phylogenetic networks
and graphs with clear separation between layout computation and rendering.
"""

from .dnetwork import plot_dnetwork
from .graphs import plot_dmgraph, plot_mmgraph
from .sdnetwork import plot_sdnetwork

# Export styles for convenience (backward compatibility)
from .dnetwork import DNetStyle as NetworkStyle, default_style

__all__ = [
    'plot_dnetwork',
    'plot_sdnetwork',
    'plot_dmgraph',
    'plot_mmgraph',
    'NetworkStyle',  # Alias for backward compatibility
    'default_style',  # Default for dnetwork (backward compatibility)
]
