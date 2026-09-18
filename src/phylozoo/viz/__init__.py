"""
Network visualization and plotting (viz).

The visualization module provides a flexible plotting system for phylogenetic
networks and graphs with support for multiple layout algorithms and customizable
styling. Use :func:`plot` for all supported types (it dispatches by object type)
and :func:`to_pretty_print` for a text drawing that needs no matplotlib.
"""

from .pretty_print import to_pretty_print
from .plot import plot

__all__ = ["plot", "to_pretty_print"]
