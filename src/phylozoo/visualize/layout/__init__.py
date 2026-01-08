"""
Layout classes for phylogenetic networks.

This module provides hierarchical layout classes for storing and managing
node positions and edge routing information for phylogenetic networks.
The layouts are renderer-agnostic and can be used with various backends.
"""

from .base import EdgeRoute, EdgeType, Layout
from .dnetwork import (
    DNetLayout,
    RectangularDNetLayout,
    compute_rectangular_dnet_layout,
    render_rectangular_dnet_layout,
)

__all__ = [
    'Layout',
    'DNetLayout',
    'RectangularDNetLayout',
    'compute_rectangular_dnet_layout',
    'render_rectangular_dnet_layout',
    'EdgeRoute',
    'EdgeType',
]

