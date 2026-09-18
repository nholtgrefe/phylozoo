"""
Network visualization and plotting (viz).

The visualization module provides a flexible plotting system for phylogenetic
networks and graphs with support for multiple layout algorithms and customizable
styling. Use :func:`plot` for all supported types (it dispatches by object type)
and :func:`to_preview_string` for a text drawing that needs no matplotlib.
"""

from .preview import to_preview_string
from .plot import plot

__all__ = ["plot", "to_preview_string"]
