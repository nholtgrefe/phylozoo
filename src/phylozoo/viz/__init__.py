"""
Network visualization and plotting (viz).

The visualization module provides a flexible plotting system for phylogenetic
networks and graphs with support for multiple layout algorithms and customizable
styling. Use :func:`plot` for all supported types (it dispatches by object type)
and :func:`to_ascii` for a text drawing that needs no matplotlib.
"""

from .ascii import to_ascii
from .plot import plot

__all__ = ["plot", "to_ascii"]
