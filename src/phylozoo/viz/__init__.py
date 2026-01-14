"""
Visualization module for PhyloZoo (viz).

This module provides a modular visualization system for phylogenetic networks
and graphs with clear separation between layout computation, edge rendering, and styling.
"""

from .dnetwork import (
    NetworkStyle,
    compute_dag_layout,
    default_style,
    get_backend,
    plot_network,
    plot_tree,
    register_backend,
)
from .graphs import plot_directed_multigraph, plot_mixed_multigraph

__all__ = [
    'plot_network',
    'plot_tree',
    'compute_dag_layout',
    'plot_directed_multigraph',
    'plot_mixed_multigraph',
    'NetworkStyle',
    'default_style',
    'get_backend',
    'register_backend',
]

