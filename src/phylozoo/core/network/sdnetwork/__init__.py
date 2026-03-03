"""
Semi-directed and mixed network module.

This module provides the SemiDirectedPhyNetwork and MixedPhyNetwork classes and
related functions for working with semi-directed phylogenetic networks. Semi-directed
networks contain both directed edges (reticulation arcs) and undirected edges (tree
edges). Internal nodes have degree >= 3, and hybrid nodes have at least one incoming
directed edge. The public API is re-exported here; the implementation is split
across the base, sd_phynetwork, features, classifications, transformations,
derivations, conversions, isomorphism, and io submodules.
"""

from .base import MixedPhyNetwork
from .sd_phynetwork import SemiDirectedPhyNetwork
from . import classifications, features, transformations, derivations, io, conversions, isomorphism

# Import io module to ensure format handlers are registered
from . import io  # noqa: F401

__all__ = [
    "SemiDirectedPhyNetwork",
    "MixedPhyNetwork",
    "classifications",
    "features",
    "transformations",
    "derivations",
    "io",
    "conversions",
    "isomorphism",
]

