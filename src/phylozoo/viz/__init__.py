"""
Network visualization and plotting (viz).

The visualization module provides a flexible plotting system for phylogenetic
networks and graphs with support for multiple layout algorithms and customizable
styling. Use :func:`plot` for all supported types. It dispatches by object type:

"""

from .plot import plot

__all__ = ['plot']
